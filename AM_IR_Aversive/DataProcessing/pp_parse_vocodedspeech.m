function [TrialData,SpoutStream,SoundStream,RateStream,Info] = pp_parse_vocodedspeech(SoundData,epocs,Info)
%
% [TrialData,SpoutStream,SoundStream,RateStream,Info] = pp_parse_vocodedspeech(SoundData,epocs,Info)
% 
% called by: ProcessPhys_SynKS
% to make TrialData table for AM aversive experiment
% 
% Currently not collecting ITIs following Warn trials
% 
% KP, 2019-04


%% StimID key for labeling in analyses

Info.stim_ID_key = { 'AsYou';  'blibBER';   'ICantBlabSuch';   'ImTheLorax';  'Please';  'Trees'  };
%   0: silence          1          2              3                  4            5         6      


%% Get Stream Data (fs = 1 kHz; 1 ms)

SpoutStream = round( resample(double(SoundData(7,:)),10000,round(Info.fs_sound*10),5) );

% Actually save downsampled RMS 
SoundStream_long = envelope(double(SoundData(2,:)),40,'rms');
SoundStream = resample(SoundStream_long,10000,round(Info.fs_sound*10),5);

RateStream = resample(double(SoundData(1,:)),10000,round(Info.fs_sound*10),5);


%%
% Onsets and offsets of trials (StimTrial TTL)
trialOnset   = round(epocs.TTyp.onset * Info.fs_sound); %samples
trialOffset  = round(epocs.TTyp.offset * Info.fs_sound); %samples
trialID      = epocs.rVID.data;
ITI_up       = 1+find(diff(SoundData(8,:))==1);

% Preallocate
Tr_Onsets  = nan(numel(trialID),1);
Tr_Offsets = nan(numel(trialID),1);
Tr_ID      = nan(numel(trialID),1);
Tr_SPL     = nan(numel(trialID),1);
Tr_LP      = nan(numel(trialID),1);

for it = 1:numel(trialID)
    
    Tr_Onsets(it)    = trialOnset(it);
    Tr_Offsets(it)   = trialOffset(it);
    
    Tr_ID(it)        = trialID(it);
    
    % Don't have to check stability
    Tr_SPL(it)       = double(SoundData(4,Tr_Onsets(it)));
    if isfield(epocs,'LPxx')
        Tr_LP(it)    = epocs.LPxx.data(it);
    else
        Tr_LP(it)    = double(SoundData(6,Tr_Onsets(it)));
    end
    
end %it


% Get timestamps of Silent period at beginning of session (for baseline)

win = round(1*Info.fs_sound);
coeffFilt   = ones(1, win)/win;
fDelay    = (length(coeffFilt)-1)/2;

Silence_Onset  = round(Info.fs_sound) + find(abs(filter(coeffFilt, 1, double(SoundData(7,:)))-1)<1e-3,1,'first')-fDelay; %when spout goes high: SoundData(7,:) + 1 second
Silence_Offset = find(SoundData(2,:),1,'first'); %when sound begins - SoundData(2,:)


%% Make TrialData
% Combine all trial data into Table (and convert samples to ms)

% Data Table [ n_trials x tags ] 
% each row is a trial, of the types listed above
% columns/tags are:  
% onset(ms)  offset(ms)  trial_ID  iti_flag  ir_flag  dB  LP  spout_flag ... behavior info ... artifact_flag

TrialData = table;
TrialData.onset   = round([Silence_Onset;  Tr_Onsets] ./Info.fs_sound.*1000); %ms
TrialData.offset  = round([Silence_Offset; Tr_Offsets]./Info.fs_sound.*1000); %ms
TrialData.trID    = [0; Tr_ID];
TrialData.SPL     = [0; Tr_SPL];
TrialData.LP      = [0; Tr_LP];

TrialData = sortrows(TrialData,'onset');



end






