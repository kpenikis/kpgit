% For each panel (pooling method) of old Fig 4, plot new d' against best SU
% d' for each stimulus.



close all

whichClass   = 'Full';

crAmt = 0.01;

%%%%%  SEPARATE ALL RESULTS BY STIMULUS


%%
% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,'eachStim');
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

% Load subpop results    
q=load(fullfile(datadir_AM,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
CR_pfr_AM = q.CR;
clear q
q=load(fullfile(datadir_AM,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
CR_Qpfr_AM = q.CR;
clear q
q=load(fullfile(datadir_AM,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
CR_dp_AM = q.CR;
clear q

% Also get CTTS
[CTTS_AM,theseCells_AM] = recallDataParams('AC','each',12);


%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitInfo_Sp = q.UnitInfo;
UnitData_Sp = q.UnitData;
clear q

% Load SU data 
q=load(fullfile(datadir_Sp,'each','CR_each.mat'));
CReach_Sp = q.CR;
clear q

% Load subpop results    
q=load(fullfile(datadir_Sp,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
CR_pfr_Sp = q.CR;
clear q
q=load(fullfile(datadir_Sp,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
CR_Qpfr_Sp = q.CR;
clear q
q=load(fullfile(datadir_Sp,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
CR_dp_Sp = q.CR;
clear q

% Also get CTTS
[CTTS_Sp,theseCells_Sp] = recallDataParams('Speech','each',12);


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
sqsmall   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];

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


%==========================================================================
% RS indices
iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);

% Significance check
% UnSig_AM = bootstrap4significance(CReach_AM(iRS_AM,:));
% UnSig_Sp = bootstrap4significance(CReach_Sp(iRS_Sp,:));

% Sort units by peakFR and dprime
% AM
[pkFR_AM,ipkFR_AM]   = rankPeakFR(CTTS_AM(iRS_AM,:,:,:));
[~,iSUdps_AM]  = sort(CReach_AM(iRS_AM,:).dprime,'descend');
% Speech
[pkFR_Sp,ipkFR_Sp]   = rankPeakFR(CTTS_Sp(iRS_Sp,:,:,:));
[~,iSUdps_Sp]  = sort(CReach_Sp(iRS_Sp,:).dprime,'descend');

% Relationship between peak FR and SU d'
figure;
plot(pkFR_AM,CReach_AM(iRS_AM(ipkFR_AM),:).dprime,'o')
[r,p]=corr(pkFR_AM,CReach_AM(iRS_AM(ipkFR_AM),:).dprime,'Type','Spearman')
hold on
plot(pkFR_Sp,CReach_Sp(iRS_Sp(ipkFR_Sp),:).dprime,'ok')
[r,p]=corr(pkFR_Sp,CReach_Sp(iRS_Sp(ipkFR_Sp),:).dprime,'Type','Spearman')



%% Get individual stimulus d'

nStim    = size(CReach_AM(1,:).Results{:},1);

% SU
dpSU_AM = nan(size(CReach_AM,1),nStim);
for ii=1:size(CReach_AM,1)    
    dpSU_AM(ii,:) = dp_from_ConfMat(mean(CReach_AM(ii,:).Results{:},3,'omitnan'),crAmt);
end
dpSU_Sp = nan(size(CReach_Sp,1),nStim);
for ii=1:size(CReach_Sp,1)    
    dpSU_Sp(ii,:) = dp_from_ConfMat(mean(CReach_Sp(ii,:).Results{:},3,'omitnan'),crAmt);
end

% Quantiles
dp_Q_AM = nan(size(CR_Qpfr_AM,1),nStim);
for ii=1:size(CR_Qpfr_AM,1)    
    dp_Q_AM(ii,:) = dp_from_ConfMat(mean(CR_Qpfr_AM(ii,:).Results{:},3,'omitnan'),crAmt);
end
dp_Q_Sp = nan(size(CR_Qpfr_Sp,1),nStim);
for ii=1:size(CR_Qpfr_Sp,1)    
    dp_Q_Sp(ii,:) = dp_from_ConfMat(mean(CR_Qpfr_Sp(ii,:).Results{:},3,'omitnan'),crAmt);
end

% Increasing pools
dp_P_AM = nan(size(CR_dp_AM,1),nStim);
for ii=1:size(CR_dp_AM,1)    
    dp_P_AM(ii,:) = dp_from_ConfMat(mean(CR_dp_AM(ii,:).Results{:},3,'omitnan'),crAmt);
end
dp_P_Sp = nan(size(CR_dp_Sp,1),nStim);
for ii=1:size(CR_dp_Sp,1)    
    dp_P_Sp(ii,:) = dp_from_ConfMat(mean(CR_dp_Sp(ii,:).Results{:},3,'omitnan'),crAmt);
end


%% For each panel (pooling method) of old Fig 4, 
%  plot new d' against best SU d' for each stimulus

ymaxval =  5.5;

% Plot
hf=figure;
set(hf,'Position',sqsmall)

% ::::::::::::::  AM  ::  Peak FR  ::  Quantiles  ::::::::::::::
subplot(2,2,3)
hold on
plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
xlabel('d'' best SU')
ylabel('d'' of pool')
box off
set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
% title('AM, quantiles, sort peak FR')

% ::::::::::::::  AM  ::  dprime  ::  SU data  ::::::::::::::
subplot(2,2,1)
hold on
plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
xlabel('d'' best SU')
ylabel('d'' of pool')
box off
set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
title('AM')

% ::::::::::::::  Speech  ::  Peak FR  ::  Quantiles  ::::::::::::::
subplot(2,2,4)
hold on
plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
xlabel('d'' best SU')
ylabel('d'' of pool')
box off
set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
% title('Speech, quantiles, sort peak FR')

% ::::::::::::::  Speech  ::  dprime  ::  SU data  ::::::::::::::
subplot(2,2,2)
hold on
plot([-0.2 ymaxval],[-0.2 ymaxval],'Color',0.32*[1 1 1])
xlabel('d'' best SU')
ylabel('d'' of pool')
box off
set(gca,'Color','none','xtick',[0 2 4],'ytick',[0 2 4])
axis square
xlim([-0.1 ymaxval])
ylim([-0.1 ymaxval])
title('Speech')


for ist = 1:8
        
    % ::::::::::::::  AM  ::  Peak FR  ::  Quantiles  ::::::::::::::
    
    subplot(2,2,3)
    hold on
    
    % Indentify SUs in each pool
    rankedSUdp = dpSU_AM(iRS_AM(ipkFR_AM),ist);
    theseSU = (CR_Qpfr_AM.iC-[zeros(size(CR_Qpfr_AM.iC,1),1)  CR_Qpfr_AM.nC-1])';
    
    for iq = 1:size(theseSU,2)
        plot( max(rankedSUdp(theseSU(2,iq):theseSU(1,iq))), dp_Q_AM(iq,ist),...
            'o','MarkerEdgeColor','none','MarkerSize',5.5,'MarkerFaceColor',colors(ist,:))
    end
    
    
    % ::::::::::::::  AM  ::  dprime  ::  SU data  ::::::::::::::
    
    subplot(2,2,1)
    hold on
    
    % Indentify SUs in each pool
    rankedSUdp = dpSU_AM(iRS_AM(iSUdps_AM),ist);
    
    N = histcounts(CR_dp_AM.iC,'BinMethod','integers');
    iCs_dp_AM = find(N>2);
    thisNC = 30;
    
    for ic = 1:numel(iCs_dp_AM)
        CRidx = find(ismember(CR_dp_AM.iC,iCs_dp_AM(ic)) & ismember(CR_dp_AM.nC,thisNC));
        plot( max(rankedSUdp(iCs_dp_AM(ic)+(1:thisNC)-1)), dp_P_AM(CRidx,ist),...
            'o','MarkerEdgeColor','none','MarkerSize',5.5,'MarkerFaceColor',colors(ist,:))
    end
    
    
    % ::::::::::::::  Speech  ::  Peak FR  ::  Quantiles  ::::::::::::::
    
    subplot(2,2,4)
    hold on
    
    % Indentify SUs in each pool
    rankedSUdp = dpSU_Sp(iRS_Sp(ipkFR_Sp),ist);
    theseSU = (CR_Qpfr_Sp.iC-[zeros(size(CR_Qpfr_Sp.iC,1),1)  CR_Qpfr_Sp.nC-1])';
    
    for iq = 1:size(theseSU,2)
        plot( max(rankedSUdp(theseSU(2,iq):theseSU(1,iq))), dp_Q_Sp(iq,ist),...
            'o','MarkerEdgeColor','none','MarkerSize',5.5,'MarkerFaceColor',colors(ist,:))
    end
    
    
    % ::::::::::::::  Speech  ::  dprime  ::  SU data  ::::::::::::::
    
    subplot(2,2,2)
    hold on
    
    % Indentify SUs in each pool
    rankedSUdp = dpSU_Sp(iRS_Sp(iSUdps_Sp),ist);
    
    N = histcounts(CR_dp_Sp.iC,'BinMethod','integers');
    iCs_dp_Sp = find(N>2);
    thisNC = 30;
    
    for ic = 1:numel(iCs_dp_Sp)
        CRidx = find(ismember(CR_dp_Sp.iC,iCs_dp_Sp(ic)) & ismember(CR_dp_Sp.nC,thisNC));
        plot( max(rankedSUdp(iCs_dp_Sp(ic)+(1:thisNC)-1)), dp_P_Sp(CRidx,ist),...
            'o','MarkerEdgeColor','none','MarkerSize',5.5,'MarkerFaceColor',colors(ist,:))
    end
    
    
end


print_eps_kp(gcf,fullfile(savedir,'BestSU_Pool_dp_eachStim'))




%% Plot pool results like old Fig 4, but one fig per stim

ymaxval = 5;

for ist = 1:8
    
    % Plot
    hf=figure;
    set(hf,'Position',sqsmall)
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Peak FR  ::  Quantiles
    
    NC_pfr_AM = mode(CR_pfr_AM.nC);
    iCRs_pfr_AM = find(CR_pfr_AM.nC==NC_pfr_AM);
    NC_dp_AM = mode(CR_dp_AM.nC);
    iCRs_dp_AM = find(CR_dp_AM.nC==NC_dp_AM);
    
    subplot(2,2,3)
    hold on
    plot( 1:length(iRS_AM),dpSU_AM(iRS_AM(ipkFR_AM),ist),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
    
    plot( (CR_Qpfr_AM.iC-[zeros(size(CR_Qpfr_AM.iC,1),1)  CR_Qpfr_AM.nC-1])', [dp_Q_AM(:,ist) dp_Q_AM(:,ist)]',...
        '-g','MarkerSize',25,'LineWidth',3)
    
    xlabel('Cell N')
    ylabel('d'' across stimuli')
    grid on
    box off
    set(gca,'Color','none','xtick',round(linspace(0,length(ipkFR_AM),6)))
    xlim([0 numel(iRS_AM)])
    ylim([-0.1 ymaxval])
    title(['AM stim ' num2str(ist) ', pkFR'])%, NC=' num2str(NC_pfr_AM)])
    
    
    % dprime  ::  SU data
    
    subplot(2,2,1)
    hold on
    plot( 1:length(iRS_AM),dpSU_AM(iRS_AM(iSUdps_AM),ist),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
    
    xlabel('Cell N')
    ylabel('d'' across stimuli')
    grid on
    box off
    set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps_AM),6)))
    xlim([0 numel(iRS_AM)])
    ylim([-0.1 ymaxval])
    title(['AM stim ' num2str(ist) ', best d''']) %, NC=' num2str(NC_dp_AM)])
    
    
    %% ~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Peak FR  ::  Quantiles
    
    NC_pfr_Sp = mode(CR_pfr_Sp.nC);
    iCRs_pfr_Sp = find(CR_pfr_Sp.nC==NC_pfr_Sp);
    NC_dp_Sp = mode(CR_dp_Sp.nC);
    iCRs_dp_Sp = find(CR_dp_Sp.nC==NC_dp_Sp);
    
    subplot(2,2,4)
    hold on
    plot( 1:length(iRS_Sp),dpSU_Sp(iRS_Sp(ipkFR_Sp),ist),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
    
    plot( (CR_Qpfr_Sp.iC-[zeros(size(CR_Qpfr_Sp.iC,1),1)  CR_Qpfr_Sp.nC-1])', [dp_Q_Sp(:,ist) dp_Q_Sp(:,ist)]',...
        '-g','MarkerSize',25,'LineWidth',3)
    
    xlabel('Cell N')
    ylabel('d'' across stimuli')
    grid on
    box off
    set(gca,'Color','none','xtick',round(linspace(0,length(ipkFR_Sp),6)))
    xlim([0 numel(iRS_Sp)])
    ylim([-0.1 ymaxval])
    title(['Speech stim ' num2str(ist) ', pkFR'])%, NC=' num2str(NC_pfr_Sp)])
    
    
    % dprime  ::  SU data
    
    subplot(2,2,2)
    hold on
    plot( 1:length(iRS_Sp),dpSU_Sp(iRS_Sp(iSUdps_Sp),ist),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
    
    xlabel('Cell N')
    ylabel('d'' across stimuli')
    grid on
    box off
    set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps_Sp),6)))
    xlim([0 numel(iRS_Sp)])
    ylim([-0.1 ymaxval])
    title(['Speech stim ' num2str(ist) ', best d''']) %, NC=' num2str(NC_dp_Sp)])
    
    
    
    %=========================
    % Add pool data (rank d')
    %=========================
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    N = histcounts(CR_dp_AM.iC,'BinMethod','integers');
    iCs_dp_AM = find(N>2);
    
    subplot(2,2,1)
    hold on
    for ic = 1:numel(iCs_dp_AM)
        CRidx = find(ismember(CR_dp_AM.iC,iCs_dp_AM(ic)));
        [ncs,incs] = sort(CR_dp_AM(CRidx,:).nC);
        
        plot(iCs_dp_AM(ic)-1+ncs,dp_P_AM(CRidx(incs),ist),'-b','LineWidth',3)
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    
    N = histcounts(CR_dp_Sp.iC,'BinMethod','integers');
    iCs_dp_Sp = find(N>4);
    
    subplot(2,2,2)
    hold on
    for ic = 1:numel(iCs_dp_Sp)
        CRidx = find(ismember(CR_dp_Sp.iC,iCs_dp_Sp(ic)));
        [ncs,incs] = sort(CR_dp_Sp(CRidx,:).nC);
        
        plot(iCs_dp_Sp(ic)-1+ncs,dp_P_Sp(CRidx(incs),ist),'-b','LineWidth',3)
    end
    
    
    print_eps_kp(gcf,fullfile(savedir,['SUvPools_all_stim' num2str(ist)]))
    
    
end %ist


%==========================================================================
%% FIG 3 
%  sort d', cumulative N spikes
%==========================================================================

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

[dps_AM,iSUdps_AM] = sort(CReach_AM(iRS_AM,:).dprime,'ascend');


idp1 = find(dps_AM<prctile(dps_AM,94),1,'last');
avgNspks = mean(sum(mean(CTTS_AM(iRS_AM,:,:,:),3,'omitnan'),2),4);
pBest_nSpk = (sum(avgNspks(iSUdps_AM((idp1+1):end)))/sum(avgNspks) + sum(avgNspks(iSUdps_AM((idp1):end)))/sum(avgNspks)) /2;
pRest_nSpk = 1-pBest_nSpk;

ymaxval = 2.5;
figure; 

yyaxis right
hold on
plot( find(~UnSig_AM(iSUdps_AM)),dps_AM(~UnSig_AM(iSUdps_AM)),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_AM(iSUdps_AM)),dps_AM(UnSig_AM(iSUdps_AM)),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])

plot([idp1 idp1]+0.5,[-0.1 ymaxval],'-','Color',0.7*[1 1 1])
ylabel('dprime')
xlim([0 numel(dps_AM)+1])
ylim([-0.1 ymaxval])

yyaxis left
plot(1:numel(dps_AM),cumsum(avgNspks(iSUdps_AM))./sum(avgNspks(iSUdps_AM)),'-','LineWidth',3)
hold on
plot([0 idp1+0.5],[pRest_nSpk pRest_nSpk],'-','Color',0.7*[1 1 1])
ylabel('cumulative proportion of spikes')
ylim([0 1])
set(gca,'Color','none','ytick',0:0.1:1,'xtick',round((0:0.2:1)*numel(dps_AM)))

xlim([0 numel(dps_AM)+1])
xlabel('Ranked cells')
title(sprintf('%0.1f, AM',pRest_nSpk*100))


print_eps_kp(gcf,fullfile(savedir,'SU_dp_propSpikes_AM'))



%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

[dps_Sp,iSUdps_Sp] = sort(CReach_Sp(iRS_Sp,:).dprime,'ascend');


idp1 = find(dps_Sp<prctile(dps_Sp,93.5),1,'last');
avgNspks = mean(sum(mean(CTTS_Sp(iRS_Sp,:,:,:),3,'omitnan'),2),4);
pBest_nSpk = (sum(avgNspks(iSUdps_Sp((idp1+1):end)))/sum(avgNspks) + sum(avgNspks(iSUdps_Sp((idp1):end)))/sum(avgNspks)) /2;
pRest_nSpk = 1-pBest_nSpk;


ymaxval = 3;
figure; 

yyaxis right
hold on
plot( find(~UnSig_Sp(iSUdps_Sp)),dps_Sp(~UnSig_Sp(iSUdps_Sp)),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_Sp(iSUdps_Sp)),dps_Sp(UnSig_Sp(iSUdps_Sp)),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])

plot([idp1 idp1]+0.5,[-0.1 ymaxval],'-','Color',0.7*[1 1 1])
ylabel('dprime')
xlim([0 numel(dps_Sp)+1])
ylim([-0.1 ymaxval])

yyaxis left
plot(1:numel(dps_Sp),cumsum(avgNspks(iSUdps_Sp))./sum(avgNspks(iSUdps_Sp)),'-','LineWidth',3)
hold on
plot([0 idp1+0.5],[pRest_nSpk pRest_nSpk],'-','Color',0.7*[1 1 1])
ylabel('cumulative proportion of spikes')
ylim([0 1])
set(gca,'Color','none','ytick',0:0.1:1,'xtick',round((0:0.2:1)*numel(dps_Sp)))

xlim([0 numel(dps_Sp)+1])
xlabel('Ranked cells')
title(sprintf('%0.1f, Speech',pRest_nSpk*100))


print_eps_kp(gcf,fullfile(savedir,'SU_dp_propSpikes_Speech'))



%% Compare ActVec vs Projection classifier results


% SU d'

% Load Proj (Full) class SU data 
q=load(fullfile(fn.figs,'ClassAM','AC','Full','each','CR_each.mat'));
CRe_Proj_AM = q.CR;
clear q

q=load(fullfile(fn.figs,'ClassSpeech','Speech','Full','each','CR_each.mat'));
CRe_Proj_Sp = q.CR;
clear q


% PLOT

hf2=figure;
set(hf2,'Position',sqsmall)

subplot(2,2,1)
plot(CRe_Proj_AM.dprime,CReach_AM.dprime,'.k')
hold on
plot([-0.2 3],[-0.2 3],'-')
axis square
set(gca,'Color','none','ytick',0:3,'xtick',0:3)
xlim([-0.2 3])
ylim([-0.2 3])

[r,p] = corr(CRe_Proj_AM.dprime,CReach_AM.dprime,'type','Spearman')
p_wsr_AM = signrank(CRe_Proj_AM.dprime,CReach_AM.dprime)
median(CReach_AM.dprime-CRe_Proj_AM.dprime)

statstext = sprintf('SINUSOIDAL\nPearson r=%0.2f, p<0.001\nsignrank p=%0.2e\nd'' ActVec-Proj median = %0.2f',r,p_wsr_AM,median(CReach_AM.dprime-CRe_Proj_AM.dprime));

subplot(2,2,3)
text(0.3,0.5,statstext)


subplot(2,2,2)
plot(CRe_Proj_Sp.dprime,CReach_Sp.dprime,'.k')
hold on
plot([-0.2 3],[-0.2 3],'-')
axis square
set(gca,'Color','none','ytick',0:3,'xtick',0:3)
xlim([-0.2 3])
ylim([-0.2 3])

[r,p] = corr(CRe_Proj_Sp.dprime,CReach_Sp.dprime,'type','Spearman')

p_wsr_Sp = signrank(CRe_Proj_Sp.dprime,CReach_Sp.dprime)
median(CReach_Sp.dprime-CRe_Proj_Sp.dprime)

statstext = sprintf('SPEECH\nPearson r=%0.2f, p<0.001\nsignrank p=%0.2e\nd'' ActVec-Proj median = %0.2f',r,p_wsr_Sp,median(CReach_Sp.dprime-CRe_Proj_Sp.dprime));

subplot(2,2,4)
text(0.3,0.5,statstext)


print_eps_kp(gcf,fullfile(savedir,'ActVec_vs_Proj_SUdp'))









