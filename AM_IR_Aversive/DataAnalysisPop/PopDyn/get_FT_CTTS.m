function [SpikesTrials,StimTrials] = get_FT_CTTS(UnitData,iUn,TrialFlags)
% called by gatherFullTrialCTTS_Speech
%
% KP, 2020-04
%

global fn trN TD Stimuli spkshift

padPreT0 = 500;
padPost  = 500;

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

SpikesTrials  = nan(1,max(TD)+padPreT0+padPost,trN,numel(Stimuli));
StimTrials    = nan(1,max(TD)+padPreT0+padPost,trN,numel(Stimuli));

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
if strcmp(TrialFlags,'sim')
    [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,unique(vertcat(Info.artifact(:).trials)),UnitData(iUn).spl,UnitData(iUn).lpn,1);
elseif strcmp(TrialFlags,'all')
    [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);
end

allStim = unique(TrialData.trID(all_TDidx))';
% if numel(allStim)<6
%     return
% end

for ist = 1:numel(Stimuli)
    
    stid = Stimuli(ist);
    
    
    %% Collect trial indices and timestamps
    
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
    
    if numel(TDidx)>trN
        keyboard
    end
    if isempty(TDidx)
        continue
    end
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
    t3 = t2 + TD(ist);
    
    
    %% Get spiking activity for each trial
    
%     figure; hold on
%     title(stid)
    for it = 1:numel(TDidx)
        
        WinBeg = t2(it)-padPreT0;
        WinEnd = t3(it)+padPost;
        
        if WinEnd > length(SoundStream)
            continue
        end
        
        sp=[]; sp = spiketimes( spiketimes>WinBeg ...
            & spiketimes<=WinEnd ) - WinBeg;
        
        SpikesTrials(1,1:length((WinBeg+1):WinEnd),it,ist) = 0;
        SpikesTrials(1,sp,it,ist) = 1;
        
        StimTrials(1,1:length((WinBeg+1):WinEnd),it,ist)   = SoundStream((WinBeg+1):WinEnd);
        
        % Plot stim to check
%         figure; hold on
%         plot(t2(it)+(0:6000)-t2(it),SoundStream(t2(it)+(0:6000)))
%         plot([t2(it) t2(it)]-t2(it),[0 0.03],'k','LineWidth',2)
%         plot([t2(it) t2(it)]+Duration-t2(it),[0 0.03],'LineWidth',2)
%         plot([s2 s2]-t2(it),[0 0.03],'k','LineWidth',2)
%         plot([s2 s2]+Duration-t2(it),[0 0.03],'LineWidth',2)
        
    end %it
    
end %ist

catch
    keyboard
end

% figure; 
% plot(permute(mean(StimTrials(1,501:1000,:,1:8),3,'omitnan'),[2 4 1 3]),...
%     'LineWidth',3)

% includethiscell = 1;


end

