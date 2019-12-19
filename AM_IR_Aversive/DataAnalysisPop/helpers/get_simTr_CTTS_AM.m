function [SpikesTrials,StimTrials,includethiscell] = get_simTr_CTTS_AM(UnitData,iUn)
% called by cumulativeSpikeCount
% called by gatherCellTimeTrialStim
%
% KP, 2019-12-16
%

global fn trN Duration Stimuli spkshift

includethiscell = 0;

padPreT0 = 500;
padPost  = 500;

if ~exist('Stimuli','var')
    Stimuli   = 1:6;
end
if ~exist('Duration','var')
    Duration  = 1000;
end

% Load data files
try
    fprintf('Loading %s sess %s...\n',UnitData(iUn).Subject,UnitData(iUn).Session)
    clear TrialData Clusters Spikes Info
    filename = sprintf( '%s_sess-%s_TrialData',UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Info',     UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    
    % Get spiketimes and shift based on calculated integration time
    if exist('Spikes','var')                                 % >>> UMS <<<
        spiketimes = unique(round(Spikes.sorted(UnitData(iUn).Channel).spiketimes(Spikes.sorted(UnitData(iUn).Channel).assigns==UnitData(iUn).Clu') * 1000 - spkshift));  %ms
    elseif exist('Clusters','var')                            % >>> KS <<<
        iClu = find([Clusters.maxChannel] == UnitData(iUn).Channel & [Clusters.clusterID] == UnitData(iUn).Clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
    end
catch
    keyboard
end


%%
try
    
SpikesTrials  = nan(1,Duration+padPreT0+padPost,trN,numel(Stimuli)+2);
StimTrials    = nan(1,Duration+padPreT0+padPost,trN,numel(Stimuli)+2);

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,unique(vertcat(Info.artifact(:).trials)),UnitData(iUn).spl,UnitData(iUn).lpn,1);
allStim = unique(TrialData.trID(all_TDidx))';

for ist = 1:numel(Stimuli)
    
    stid = Stimuli(ist);
    
    Silence = false;
    if stid==9
        Silence = true;
    elseif ~ismember(stid,allStim)
        continue
    end
    
    
    %% Collect trial indices and timestamps
    
    if ~Silence                                          % SOUND TRIALS
        
        if stid==3 || stid==6
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
            % Find Pdc trials that follow same rate during ITI
            TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
            
            TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
            TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        else
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
            TDidx_iti = [];
        end
        
        % Get timestamps of onsets and offsets
        clear t2 t3 t_win
        t2 = TrialData.onset(TDidx);
        
        % Add ITI trials (shortened to match duration)
        if ~isempty(TDidx_iti) && numel(t2)<trN
            t2 = [t2; TrialData.onset(TDidx_iti)];
            TDidx = [TDidx; TDidx_iti];
        end
        
    else                                               % SILENCE TRIALS
        keyboard
        clear t2 t3 t_win
        SilPd = [TrialData.onset(1) TrialData.offset(1)];
        
        t2 = TrialData.onset(1) : Duration : (TrialData.offset(1)-mod(diff(SilPd),Duration)-Duration);
        TDidx = 1:length(t2);
        Ntrials(stid) = length(t2);
        
    end
    
%     if numel(t2)<trN
%         sprintf('not enough trials for this stim/cell')
%         includethiscell = 0;
%         return
%     end
    
    kt     = 1:length(t2); %randperm(length(t2),trN);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    if Silence && t3(end)>TrialData.offset(1)
        keyboard
    end
    
    
    %% Get spiking activity for each trial
    
    for it = 1:min(numel(TDidx),trN)
        
        istid = stid;
        if stid==8
            istid = 9; %advance one slot for DB
        end
        
        WinBeg = t2(it)-padPreT0;
        WinEnd = t3(it)+padPost;
        
        sp=[]; sp = spiketimes( spiketimes>WinBeg ...
            & spiketimes<=WinEnd ) - WinBeg;
        
        SpikesTrials(1,:,it,istid)  = 0;
        SpikesTrials(1,sp,it,istid) = 1;
        
        StimTrials(1,:,it,istid)    = SoundStream((WinBeg+1):WinEnd);
        
        
        % Next segment start
        if stid>6
            istid = istid+1;
            
            s2 = Phase0(:,find(Phase0(1,:)>=(t2(it)+500),1,'first'));
            
            WinBeg = s2(1)-padPreT0;
            WinEnd = s2(1)+Duration+padPost;
            
            sp=[]; sp = spiketimes( spiketimes>WinBeg ...
                & spiketimes<=WinEnd ) - WinBeg;
            
            SpikesTrials(1,:,it,istid)  = 0;
            SpikesTrials(1,sp,it,istid) = 1;
            
            StimTrials(1,:,it,istid)    = SoundStream((WinBeg+1):WinEnd);
        end
        
%         % Next segment start
%         if stid==8
%             istid = istid+1;
%             s3 = Phase0(:,find(Phase0(1,:)>(s2(1)+500),1,'first'));
%             
%             sp=[]; sp = spiketimes( spiketimes>s3(1) ...
%                 & spiketimes<=(s3(1)+Duration) ) - s3(1);
%             
%             SpikesTrials(1,:,it,istid)  = 0;
%             SpikesTrials(1,sp,it,istid) = 1;
%         end
%         figure; 
%         plot(t2(it)+(0:2000),SoundStream(t2(it)+(0:2000)))
%         hold on
%         plot([t2(it) t2(it)],[0 1.5],'k','LineWidth',2)
%         plot([s2(1) s2(1)],[0 1.5])
%         plot([s3(1) s3(1)],[0 1.5])
%         plot([s3(1) s3(1)]+500,[0 1.5])
%         plot([t2(it) t2(it)]+1937,[0 1.5],'k','LineWidth',2)
        
    end %it
    
end %ist

catch
    keyboard
end

includethiscell = 1;


end

