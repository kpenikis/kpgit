function StimDPs_RandPop
% StimDPs_RandPop
%
% Plot d' of each random pool against the best SU member, for each stimulus.
%
% KP, 2020-04
%

% close all

whichClass   = 'Full';
whichCells   = 'Rand_RS';

crAmt = 0.01;


%%
% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,'eachStim',whichCells);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass);
q = load(fullfile(fn.processed,'Units'));
UnitInfo_AM = q.UnitInfo;
UnitData_AM = q.UnitData;
clear q

% Load SU data 
q=load(fullfile(datadir_AM,'each','CR_each.mat'));
CReach_AM = q.CR;
clear q

% Load rand subpop results    
q=load(fullfile(datadir_AM,whichCells,sprintf('CR_%s_%s.mat',whichClass,whichCells)));
CRpop_Cat_AM = q.CR;
clear q
q=load(fullfile(datadir_AM,whichCells,'Sum',sprintf('CR_%s_%s.mat',whichClass,whichCells)));
CRpop_Sum_AM = q.CR;
clear q

% Also get CTTS
[CTTS_AM,theseCells_AM] = recallDataParams('AC','each',12);


%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

% datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
% q = load(fullfile(fn.processed,'UnitsVS'));
% UnitInfo_Sp = q.UnitInfo;
% UnitData_Sp = q.UnitData;
% clear q
% 


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
sqsmall   = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

% Stimulus colors
colors = [ 150 150 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];


ymaxval =  5;

%==========================================================================
% RS indices
iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
% iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);


%% Get individual stimulus d'

nStim    = size(CReach_AM(1,:).Results{:},1);

% SU
dpSU_AM = nan(size(CReach_AM,1),nStim);
for ii=1:size(CReach_AM,1)    
    dpSU_AM(ii,:) = dp_from_ConfMat(mean(CReach_AM(ii,:).Results{:},3,'omitnan'),crAmt);
end
% dpSU_Sp = nan(size(CReach_Sp,1),nStim);
% for ii=1:size(CReach_Sp,1)    
%     dpSU_Sp(ii,:) = dp_from_ConfMat(mean(CReach_Sp(ii,:).Results{:},3,'omitnan'),crAmt);
% end

% Random populations
dpPop_Cat_AM     = nan(size(CRpop_Cat_AM,1),nStim);
iCRs_Pop_Cat_AM  = nan(size(CRpop_Cat_AM,1),CRpop_Cat_AM.nC(1));
dpPop_BestSU_AM  = nan(size(CRpop_Cat_AM,1),1);
dpPop_medianSU_AM  = nan(size(CRpop_Cat_AM,1),1);
dpPop_meanSU_AM  = nan(size(CRpop_Cat_AM,1),1);
for ii=1:size(CRpop_Cat_AM,1)    
    dpPop_Cat_AM(ii,:)    = dp_from_ConfMat(mean(CRpop_Cat_AM(ii,:).Results{:},3,'omitnan'),crAmt);
    iCRs_Pop_Cat_AM(ii,:) = CRpop_Cat_AM(ii,:).CRids{:}';
    dpPop_BestSU_AM(ii)   = max(CReach_AM.dprime(iCRs_Pop_Cat_AM(ii,:)));
    dpPop_medianSU_AM(ii) = median(CReach_AM.dprime(iCRs_Pop_Cat_AM(ii,:)));
    dpPop_meanSU_AM(ii)   = mean(CReach_AM.dprime(iCRs_Pop_Cat_AM(ii,:)));
end

dpPop_Sum_AM     = nan(size(CRpop_Sum_AM,1),nStim);
iCRs_Pop_Sum_AM  = nan(size(CRpop_Sum_AM,1),CRpop_Sum_AM.nC(1));
for ii=1:size(CRpop_Sum_AM,1)    
    dpPop_Sum_AM(ii,:)  = dp_from_ConfMat(mean(CRpop_Sum_AM(ii,:).Results{:},3,'omitnan'),crAmt);
    iCRs_Pop_Sum_AM(ii,:)   = CRpop_Sum_AM(ii,:).CRids{:}';
end


% Get max SU dp for each population for each stimulus
maxSUdps = nan(size(iCRs_Pop_Cat_AM,1),8);
for ist = 1:8
    for ip = 1:size(iCRs_Pop_Cat_AM,1)
        maxSUdps(ip,ist) = max(dpSU_AM(iCRs_Pop_Cat_AM(ip,:),ist));
    end
end



%% PLOT

% ::::::::::::::::::::::::::   AM   ::::::::::::::::::::::::::

% Avg d'
hfavg=figure;
set(hfavg,'Position',sqsmall)

plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
hold on
plot(dpPop_BestSU_AM,CRpop_Cat_AM.dprime,'.k')
% plot(dpPop_BestSU_AM-dpPop_meanSU_AM,CRpop_Cat_AM.dprime-dpPop_BestSU_AM,'.k')
% plot(dpPop_meanSU_AM-dpPop_medianSU_AM,CRpop_Cat_AM.dprime-dpPop_BestSU_AM,'.k')

set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
box off
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
xlabel('d'' best SU')
ylabel('d'' of random pool')
title('Avg d'' for task')


% All stimuli together
hfallam=figure;
set(hfallam,'Position',sqsmall)

plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
hold on
plot(maxSUdps(:),dpPop_Cat_AM(:),'.k')

set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
box off
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
xlabel('d'' best SU')
ylabel('d'' of random pool')
title('All stimuli')

print_eps_kp(hfallam,fullfile(savedir,'AC_All_BestvPool'))



% Each stimulus separately

hfam=figure;
set(hfam,'Position',fullscreen)

for ist = 1:8
    
    subplot(2,4,ist)
    plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
    hold on
    
    plot(maxSUdps(:,ist),dpPop_Cat_AM(:,ist),'.k')
    
    set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
    axis square
    box off
    xlim([-0.1 ymaxval])
    ylim([-0.1 ymaxval])
    xlabel('d'' best SU')
    ylabel('d'' of random pool')
    title(ist)
    
end

print_eps_kp(hfam,fullfile(savedir,'AC_Each_BestvPool'))



%% Distribution of d' for each stimulus 

hf=figure;
set(hf,'Position',fullscreen)

for ist = 1:8
    
    plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
    
    plotSpread(dpPop_Cat_AM(:,ist),'distributionIdx',ist*ones(size(dpPop_Cat_AM(:,ist))),'showMM',3)
    
end


%% Compare Cat and Sum distributions

figure;
plotSpread(CRpop_Cat_AM.dprime,'distributionIdx',1*ones(size(CRpop_Cat_AM.dprime)),'showMM',3)
hold on
plotSpread(CRpop_Sum_AM.dprime,'distributionIdx',2*ones(size(CRpop_Sum_AM.dprime)),'showMM',3)






