function gatherCellTimeTrialStim_Speech
% gatherCellTimeTrialStim_Speech
%
%   Gathering speech data for classifier.
% 
% KP, 2019-12-16
%

global fn trN Duration Stimuli spkshift

trN      = 50;
Duration = 500;
Stimuli  = 2:6;

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%%

Cell_Time_Trial_Stim = []; 
Env_Time_Trial_Stim  = []; 
Un_Indices = [];

for iUn = 1:numel(UnitData)
    
    [SpikesTrials,StimTrials,includethiscell] = get_simTr_CTTS_Speech(UnitData,iUn);
    
    if includethiscell
        Cell_Time_Trial_Stim = [Cell_Time_Trial_Stim; SpikesTrials];
        Env_Time_Trial_Stim  = [Env_Time_Trial_Stim; StimTrials];
        Un_Indices           = [Un_Indices; iUn];
    end
    
end

savedir = fullfile(fn.figs,'ClassSpeech','RawData');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

save(fullfile(savedir,'CTTS_Speech_sim'),'Cell_Time_Trial_Stim','Env_Time_Trial_Stim','Un_Indices','-v7.3')



end