
% Load data files
try
% if (iUn>1 && ~( strcmp(UnitData(iUn).Subject,UnitData(iUn-1).Subject) && strcmp(UnitData(iUn).Session,UnitData(iUn-1).Session) )) || iUn==1
    fprintf('Loading %s sess %s...\n',UnitData(iUn).Subject,UnitData(iUn).Session)
    clear TrialData Clusters Spikes Info
    filename = sprintf( '%s_sess-%s_TrialData',UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Info',     UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
% end


% Get spiketimes and shift based on calculated integration time
                                                               % >>> KS <<<
iClu = find([Clusters.maxChannel] == UnitData(iUn).Channel & [Clusters.clusterID] == UnitData(iUn).Clu);
spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 )');

% Calculate zFR
bs_hist = 1;
bs_smth = 20;
keyboard

[Stream_FRsmooth,Stream_zscore] = convertSpiketimesToFR(spiketimes,...
    TrialData.offset(end)+100,TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');


catch
    keyboard
end

iu_FRvec   = nan(1,Duration,numel(Stimuli));
iu_zFRvec  = nan(1,Duration,numel(Stimuli));
RMS        = nan(numel(Stimuli),Duration);

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);
allStim = unique(TrialData.trID(all_TDidx))';

for ist = 1:numel(Stimuli)
    
    stid = Stimuli(ist);
    
    
    %% Collect trial indices and timestamps
    
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
    
    kt     = randperm(length(t2),min(Ntrials(Ntrials>0)));
%   kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration-1;
    
    
    %% Get spiking activity for each trial
    
    FRtr  = nan(numel(TDidx),Duration);
    zFRtr = nan(numel(TDidx),Duration);
    RMStr = nan(numel(TDidx),Duration);
    
    for it = 1:numel(TDidx)
        
        FRtr(it,:)  = Stream_FRsmooth(t2(it):t3(it));
        zFRtr(it,:) = Stream_zscore(t2(it):t3(it));
        RMStr(it,:) = SoundStream(t2(it):t3(it));
        
    end %it
    
    % Save  data
    iu_FRvec(1,:,ist)  = mean(FRtr,1);
    iu_zFRvec(1,:,ist) = mean(zFRtr,1);
    RMS(ist,:)         = mean(RMStr,1);
    
    
end %ist
