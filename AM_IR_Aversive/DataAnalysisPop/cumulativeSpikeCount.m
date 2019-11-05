function cumulativeSpikeCount


global fn trN Duration Stimuli spkshift

trN      = 20;
Duration = 1000;
Stimuli  = 1:6;

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

savedir = fullfile(fn.figs,'CmSpCount');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%

Cell_Time_Trial_Stim = []; 
Un_Indices = [];

for iUn = 1:numel(UnitData)
    
    [SpikesTrials,includethiscell] = get_rand_tr_spike_trains(UnitData,iUn);
    
    if includethiscell
        Cell_Time_Trial_Stim = [Cell_Time_Trial_Stim; SpikesTrials];
        Un_Indices           = [Un_Indices; iUn];
    end
    
end

if any(any(any(any(isnan(Cell_Time_Trial_Stim)))))
    keyboard
end

save(fullfile(savedir,['Cell_Time_Trial_Stim_' num2str(trN) 'trs']),'Cell_Time_Trial_Stim','Un_Indices')


%%


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

% Set colors
colors = [   150 150 150;...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;
colors = [ colors; ...
    [37  84 156]./255 ;...
    [19 125 124]./255 ;...
    [ 0   0   0]./255];


hf=figure;
hold on
set(hf,'Position',halfscreen)

for ist=1:numel(Stimuli)
    
    CountsTrials = cumsum(sum(Cell_Time_Trial_Stim(:,:,:,ist),1),2);
    CountsTrials = permute(CountsTrials,[2,3,1]);
    
    Count_mu = mean(CountsTrials,2);
    Count_SD = std(CountsTrials,0,2);
    patch_border = [Count_mu; flipud(Count_mu)] + [Count_SD; -1*flipud(Count_SD)];
    
    patch([1:1000 1000:-1:1], patch_border', colors(ist,:),'EdgeColor',colors(ist,:),'FaceAlpha',1)
    
end

axis square
xlim([0 500])
xlabel('Time')
ylabel('Number of spikes')
title(['Cumulative spike count across population, N=' num2str(size(Cell_Time_Trial_Stim,1))])



% Save figure
print_eps_kp(hf,fullfile(savedir,sprintf('CSC_%iuns_%itrs',size(Cell_Time_Trial_Stim,1),trN)))



end