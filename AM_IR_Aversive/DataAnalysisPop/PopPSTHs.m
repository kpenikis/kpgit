function PopPSTHs
% PopPSTHs
%  
%  Raw and smoothed FR across all cells for each stim.
%
%  KP, 2019-11
%

%% Load data

close all
fn = set_paths_directories('','',1);

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
quartscreen  = [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];

% Save dir
savedir = fullfile(fn.figs,'PopPSTH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% 

% Get data
Stimuli   = 1:8;
Duration  = 1000;
skipOnset = 0;
smoothBin = 10;

All_PSTHs = nan(numel(UnitData),Duration,numel(Stimuli));
StimPlots = nan(numel(UnitData),Duration,numel(Stimuli));

for iUn = 1:numel(UnitData)
    get_psth_posthoc
    All_PSTHs(iUn,:,:) = permute(Unit_PSTHs,[2 1 3]);
    StimPlots(iUn,:,:) = permute(Unit_STIM, [2 1 3]);
end



% Plot PSTH for each stimulus
for ist = 1:numel(Stimuli)
    
    figure;
    set(gcf,'Position',quartscreen)
    
    subplot(4,1,1)
    plot(mean(StimPlots(:,:,ist),1,'omitnan'),'k','LineWidth',3)
    set(gca,'Color','none','xtick',[],'ytick',[])
    
    subplot(4,1,2:4)
    plot(sum(All_PSTHs(:,:,ist),1,'omitnan'),'k','LineWidth',3)
    ylim([0 3000])
    set(gca,'Color','none','tickdir','out')
    
    suptitle(Info.stim_ID_key{Stimuli(ist)})
    
    savename = sprintf('PopPSTH_%s',Info.stim_ID_key{Stimuli(ist)});
    print_eps_kp(gcf,fullfile(savedir,savename))
    
end

end




