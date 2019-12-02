
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

[Stream_FRsmooth,Stream_zscore] = convertSpiketimesToFR(spiketimes,...
    TrialData.offset(end)+100,TrialData.onset(1),TrialData.offset(1),bin_smooth,bin_smooth,'silence');


catch
    keyboard
end


Stimuli = [1 5 6 2];

iu_FRvec   = nan(1,Duration,numel(Stimuli));
iu_zFRvec  = nan(1,Duration,numel(Stimuli));
RMS        = nan(numel(Stimuli),Duration);

thisRaster = nan(1,Duration,numel(Stimuli),10);

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);
allStim = unique(TrialData.trID(all_TDidx))';

for ist = 1:numel(Stimuli)
    
    stid = Stimuli(ist);
    
    SegTemplate = getfield(k,kfns{ist});
    
    
    %% Collect trial indices and timestamps
    
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
    
%     kt     = randperm(length(t2),min(Ntrials(Ntrials>0)));
    kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    t3 = TrialData.offset(TDidx);
    
    nReps = round(mean(t3-t2)./length(SegTemplate));
    
    tEdgS = repmat([0 length(SegTemplate)-1],nReps,1) + ([1:nReps]'-1)*(length(SegTemplate)-1);
    
    
    %% Get spiking activity for each trial
    
    FRtr  = nan(numel(t2)*nReps,length(SegTemplate));
    zFRtr = nan(numel(t2)*nReps,length(SegTemplate));
    RMStr = nan(numel(t2)*nReps,length(SegTemplate));
    Rastr = zeros(numel(t2)*nReps,length(SegTemplate));
    
    iseg = 0;
    for it = 1:numel(t2)
        
        for irep = 1:nReps
            iseg = iseg+1;
            
            FRtr(iseg,:)  = Stream_FRsmooth( (t2(it) + tEdgS(irep,1)) : (t2(it) + tEdgS(irep,2)) );
            zFRtr(iseg,:) = Stream_zscore( (t2(it) + tEdgS(irep,1)) : (t2(it) + tEdgS(irep,2)) );
            RMStr(iseg,:) = SoundStream( (t2(it) + tEdgS(irep,1)) : (t2(it) + tEdgS(irep,2)) );
            
            sp = spiketimes( spiketimes>(t2(it) + tEdgS(irep,1)) & spiketimes<=(t2(it) + tEdgS(irep,2)) )  - (t2(it) + tEdgS(irep,1));
            Rastr(iseg,sp) = 1; 
            
        end
    end %it
    
    % Save  data
    iu_FRvec(1,1:length(SegTemplate),ist)  = mean(FRtr,1);
    iu_zFRvec(1,1:length(SegTemplate),ist) = mean(zFRtr,1);
    RMS(ist,1:length(SegTemplate))         = mean(RMStr,1);
    
    % Save 10 random trials of spiking data
    ridx = randperm(size(Rastr,1),10);
    thisRaster(1,1:length(SegTemplate),ist,:) = Rastr(ridx,:)';
%     sp_trs(iUn,1:size(thisRaster,2),stid-1,:) = thisRaster(ridx,:)';
%     thisRaster
    
end %ist



