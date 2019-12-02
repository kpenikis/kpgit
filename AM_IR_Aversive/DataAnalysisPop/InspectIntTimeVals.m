function InspectIntTimeVals
% 
% InspectIntTimeVals
% 


%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

histedges = 0:60;


%%

hf=figure;
histogram([UnitData([UnitData.IntTime_spk]>0).IntTime_spk],histedges,'FaceColor','k')
hold on
plot(mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]),0,'.b','MarkerSize',30)
plot(median([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]),0,'.g','MarkerSize',30)
plot(geomean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]),0,'.r','MarkerSize',30)
plot(harmmean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]),0,'.m','MarkerSize',30)
plot(prctile([UnitData([UnitData.IntTime_spk]>0).IntTime_spk],25),0,'.c','MarkerSize',30)


print_eps_kp(hf,fullfile(fn.figs,'IntTimeDistribution'))


end