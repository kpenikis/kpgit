% called by PopPSTHs

if ~exist('Stimuli','var')
    Stimuli   = [9 2:6];
end
if ~exist('Duration','var')
    Duration = 750;
end
if ~exist('skipOnset','var')
    skipOnset = 250;
end


% Load data files
try
if (iUn>1 && ~( strcmp(UnitData(iUn).Subject,UnitData(iUn-1).Subject) && strcmp(UnitData(iUn).Session,UnitData(iUn-1).Session) )) || iUn==1
    fprintf('Loading %s sess %s...\n',UnitData(iUn).Subject,UnitData(iUn).Session)
    clear TrialData Clusters Spikes Info
    filename = sprintf( '%s_sess-%s_TrialData',UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Info',     UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
end

% if UnitData(iUn).IntTime_spk==0
%     return
% else
%     spkshift = UnitData(iUn).IntTime_spk;
% end

% Get spiketimes and shift based on calculated integration time
if exist('Spikes','var')                                 % >>> UMS <<<
    spiketimes = unique(round(Spikes.sorted(UnitData(iUn).Channel).spiketimes(Spikes.sorted(UnitData(iUn).Channel).assigns==UnitData(iUn).Clu') * 1000 - spkshift));  %ms
elseif exist('Clusters','var')                            % >>> KS <<<
    iClu = find([Clusters.maxChannel] == UnitData(iUn).Channel & [Clusters.clusterID] == UnitData(iUn).Clu);
    spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
end
catch
    keyboard
    return
end



FRtrials      = nan(50,numel(Stimuli));
FFstim        = nan(1,numel(Stimuli));
StimSpikeData = nan(numel(Stimuli),2);

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);
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
        t2 = TrialData.onset(TDidx) + skipOnset;
        
        % Add ITI trials (shortened to match duration)
        if ~isempty(TDidx_iti)
            t2 = [t2; TrialData.onset(TDidx_iti)];
            TDidx = [TDidx; TDidx_iti];
        end
        
    else                                               % SILENCE TRIALS
        clear t2 t3 t_win
        SilPd = [TrialData.onset(1) TrialData.offset(1)];
        
        t2 = TrialData.onset(1) : Duration : (TrialData.offset(1)-mod(diff(SilPd),Duration)-Duration);
        TDidx = 1:length(t2);
        Ntrials(stid) = length(t2);
        
    end
    
    if exist('trMax','var')
        trLim = min([trMax Ntrials(Ntrials>0)]);
    else
        trLim = min(Ntrials(Ntrials>0));
    end
    kt     = randperm(length(t2),trLim);
%     kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    if Silence && t3(end)>TrialData.offset(1)
        keyboard
    end
    
    
    %% Get spiking activity for each trial
    
    nspks = nan( numel(TDidx), 1 );
    for it = 1:numel(TDidx)
        
        sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
            & spiketimes<=t3(it) ) - t2(it) - 1;
        nspks(it)   = numel(sp);
        
    end %it
    
    % Save  data
    FRtrials(1:length(nspks),ist) = nspks./(Duration/1000);
    StimSpikeData(ist,:)          = [mean(nspks) var(nspks)];
    
    
end %ist

