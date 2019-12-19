function [PSTH_cell,STIM_cell] = getdata4TRF(UnitData,iUn)
% getdata4TRF
%
% Get PSTH and STIM cells from AM or VS data
%  PSTH_cell{st}(trials,time)
%  STIM_cell{st}(trials,time)
%
%

global fn spkshift smth_win exclOnset AM_durs VS_durs


% Get unit info
subject     = UnitData(iUn).Subject;
session     = UnitData(iUn).Session;
channel     = UnitData(iUn).Channel(1);
clu         = UnitData(iUn).Clu(1);

% Get sound parameters
dBSPL       = UnitData(iUn).spl;
LP          = UnitData(iUn).lpn;


% Load data files
fprintf('Loading %s sess %s ch %i clu %i...\n',subject,session,channel,clu)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename)); 
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));


% Get spiketimes and shift based on calculated integration time
if exist('Spikes','var')                                 % >>> UMS <<<
    
    spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift));  %ms
    %         spiketimes = unique(round(  (Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift) / 1000 * FSfs )); %samples
    
elseif exist('Clusters','var')                            % >>> KS <<<
    
    iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
    spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
    
end


% If shuffling spiketimes
%     spiketimes = sort(randi( length(SpoutStream), size(spiketimes) ));


%%

[Stream_FRsmooth,Stream_zscore,Stream_spikes,~] = convertSpiketimesToFR(spiketimes,...
    length(SoundStream),TrialData.onset(1),TrialData.offset(1),'exp',smth_win,'silence');

% Find all stimuli presented with these parameters, and trials without 
% diruptive artifact & while the animal was drinking
[all_TDidx,Ntrials,~,allStim] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP,0);
Ntrials(1) = nan; %so Warn stimulus doesn't interfere

% Set stim vars depending on what type of session
if strcmp(session(end-1:end),'AM')
%     theseStIDs = intersect(2:8,allStim');
    theseStIDs = 2:8;
    Durations = AM_durs;
elseif strcmp(session(end-1:end),'VS')
    theseStIDs = 1:6;
    Durations = VS_durs;
end

STIM_cell = cell(numel(theseStIDs),1);
PSTH_cell = cell(numel(theseStIDs),1);

% N trials of Stim and PSTH per group
% Pdc   --> Irr 1
% Irr 2 --> Irr 1  %compare to within group prediction, but need different trials

% Do need to check for number of trials, but do this later?
% Differs for Pdc-Irr vs AM-VS comparisons

% nTrGrp = min(Ntrials(theseStIDs)); %floor(min(Ntrials(2:end))/2);

% if any(Ntrials(7:end) < nTrGrp*2)
%     thisIR = allStim(Ntrials(7:end) < nTrGrp*2) + 6;
%     all_TDidx(TrialData.trID(all_TDidx)==thisIR)  = [];
%     Ntrials(allStim==thisIR) = [];
%     allStim(allStim==thisIR) = [];
% end


%%

for stid = theseStIDs
    
    ist = find(stid==theseStIDs);
    
    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==stid);
    
    if numel(st_TDidx_ALL)<8
        continue
    end
    
    %%%  Skip ITI stimuli
    if strcmp(session(end-1:end),'AM')
        ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
        TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(1));
    else
        TDidx = st_TDidx_ALL;
    end
    
    % Get timestamps of onsets and offsets
    clear t2 t3 Duration t_win
    t2 = TrialData.onset(TDidx);
    Duration = Durations(stid); %mode(diff([t2 t3],1,2));
    t3 = t2 + Duration;
    if t3(end)>length(Stream_FRsmooth)
        t3 = t3(1:end-1);
        t2 = t2(1:end-1);
    end
    
    if exclOnset
        t2 = t2+150;
    end
    
    % Collect responses for each trial
    stim_i   = nan(  numel(t2), Duration);
    psth_rw  = zeros(numel(t2), Duration);
    psth_sm  = nan(  numel(t2), Duration);
    
    for it = 1:numel(t2)
        
        stim_i(it,:) = ...
            SoundStream(1, t2(it) : (t3(it)-1) ); 
        
        sp=[]; sp = spiketimes( spiketimes>t2(it) ...
            & spiketimes<=t3(it) ) - t2(it)  ;
        
        psth_rw(it,sp) = 1000;  % raw, compare to smoothed
        psth_sm(it,:) = Stream_FRsmooth(1, t2(it) : (t3(it)-1) );
        
    end %it
    
    STIM_cell{ist} = exp( stim_i );
    PSTH_cell{ist} = psth_sm;
    
end %istim


% Normalize envelope across stimuli
stimMin = min(cell2mat(cellfun(@(x) min(min(x)),STIM_cell,'UniformOutput',false)));
foo = cellfun(@(x) x-stimMin, STIM_cell,'UniformOutput',false);
%     STIM = STIM/max(max(STIM));
STIM_cell = foo; clear foo


% OPTIONAL: Subtract spontaneous FR from PSTH
%     foo = cellfun(@(x) max(x-UnitData(iUn).BaseFR,0), PSTH_cell, 'UniformOutput',false);
%     PSTH_cell = foo; clear foo



end