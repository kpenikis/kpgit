function gatherCTTSContext_Speech
% gatherCellTimeTrialStim_Speech
%
%   Gathering speech data for classifier.
% 
% KP, 2020-01-22
%

global fn trN Duration Stimuli spkshift k kfns SegDurs

trN      = 350;


% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


% Load RepSpeech segment templates
k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
kfns = fieldnames(k);
SegDurs  = structfun(@length,k);
Duration = max(structfun(@length,k));

% 4 segments
Stimuli = [1 5 6 2];

% 2 (or 3) contexts



%%

Cell_Time_Trial_Stim = []; 
Env_Time_Trial_Stim  = []; 
Un_Indices = [];

for iUn = 1:numel(UnitData)
    
    [SpikesTrials,StimTrials,includethiscell] = get_simTr_CTTScontext_Speech(UnitData,iUn);
    
    if includethiscell
        Cell_Time_Trial_Stim = [Cell_Time_Trial_Stim; SpikesTrials];
        Env_Time_Trial_Stim  = [Env_Time_Trial_Stim;    StimTrials];
        Un_Indices           = [Un_Indices; iUn];
    end
    
end

savedir = fullfile(fn.figs,'ClassContext','Speech','RawData');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

save(fullfile(savedir,'CTTSC_Speech_sim'),'Cell_Time_Trial_Stim','Env_Time_Trial_Stim','Un_Indices','-v7.3')



end