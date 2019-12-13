function gatherCellTimeTrialStim


global fn trN Duration Stimuli spkshift

trN      = 100;
Duration = 500;
Stimuli  = 1:8;

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%%

Cell_Time_Trial_Stim = []; 
Un_Indices = [];

for iUn = 1:numel(UnitData)
    
    [SpikesTrials,includethiscell] = get_rand_tr_spike_trains_allTrs(UnitData,iUn);
    
    if includethiscell
        Cell_Time_Trial_Stim = [Cell_Time_Trial_Stim; SpikesTrials];
        Un_Indices           = [Un_Indices; iUn];
    end
    
end

savedir = fullfile(fn.figs,'StimClass');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

save(fullfile(savedir,'Cell_Time_Trial_Stim'),'Cell_Time_Trial_Stim','Un_Indices','-v7.3')



end