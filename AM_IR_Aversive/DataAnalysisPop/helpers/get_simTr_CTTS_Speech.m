function [SpikesTrials,StimTrials,includethiscell] = get_simTr_CTTS_Speech(UnitData,iUn)
% called by gatherCellTimeTrialStim_Speech
%
% KP, 2019-12-16
%

global fn trN Duration Stimuli spkshift

includethiscell = 0;

spkshift = 0;

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
% [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,unique(vertcat(Info.artifact(:).trials)),UnitData(iUn).spl,UnitData(iUn).lpn,1);
[all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);

allStim = unique(TrialData.trID(all_TDidx))';
% if numel(allStim)<6
%     return
% end

for ist = 1:numel(Stimuli)
    
    stid = Stimuli(ist);
    
    
    %% Collect trial indices and timestamps
    
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
        
%     if numel(t2)<11
% %         keyboard
%         sprintf('not enough trials for this stim/cell')
%         includethiscell = 0;
%         return
%     end
    
    kt     = 1:length(t2);
    t2     = t2(kt);
    if stid==3
        t2 = t2+30;
    end
    if stid==4
        t2 = t2+60;
    end
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    
    %% Get spiking activity for each trial
%     figure; hold on
%     title(stid)
    for it = 1:min(numel(TDidx),trN)
        
        istid = stid;
        if stid==4
            istid = 5; %advance one slot for DB
        end
        if stid>4
            istid = stid+2; %advance one slot for DB
        end
        
        
        WinBeg = t2(it)-padPreT0;
        WinEnd = t3(it)+padPost;
        
        sp=[]; sp = spiketimes( spiketimes>WinBeg ...
            & spiketimes<=WinEnd ) - WinBeg;
        
        SpikesTrials(1,:,it,istid)  = 0;
        SpikesTrials(1,sp,it,istid) = 1;
        
        StimTrials(1,:,it,istid) = SoundStream((WinBeg+1):WinEnd);
        
        
        % Add another segment from sentence STIMULUS 3
        if stid==3
            istid = istid+1;
            
            s2 = t2(it) + 650;
            
            WinBeg = s2-padPreT0;
            WinEnd = s2+Duration+padPost;
            
            sp=[]; sp = spiketimes( spiketimes>WinBeg ...
                & spiketimes<=WinEnd ) - WinBeg;
            
            SpikesTrials(1,:,it,istid)  = 0;
            SpikesTrials(1,sp,it,istid) = 1;
            
            StimTrials(1,:,it,istid) = SoundStream((WinBeg+1):WinEnd);
        end
        
        
        % Add another segment from sentence STIMULUS 4
        if stid==4 
            istid = istid+1;
            
            s2 = t2(it) + 2810;
            
            WinBeg = s2-padPreT0;
            WinEnd = s2+Duration+padPost;
            
            sp=[]; sp = spiketimes( spiketimes>WinBeg ...
                & spiketimes<=WinEnd ) - WinBeg;
            
            SpikesTrials(1,:,it,istid)  = 0;
            SpikesTrials(1,sp,it,istid) = 1;
            
            StimTrials(1,:,it,istid) = SoundStream((WinBeg+1):WinEnd);
        end
        
        % Plot stim to check
%         figure; hold on
%         plot(t2(it)+(0:2000)-t2(it),SoundStream(t2(it)+(0:2000)))
%         plot([t2(it) t2(it)]-t2(it),[0 0.03],'k','LineWidth',2)
%         plot([t2(it) t2(it)]+Duration-t2(it),[0 0.03],'LineWidth',2)
%         plot([s2 s2]-t2(it),[0 0.03],'k','LineWidth',2)
%         plot([s2 s2]+Duration-t2(it),[0 0.03],'LineWidth',2)
        
    end %it
    
end %ist

catch
    keyboard
end

includethiscell = 1;


end

