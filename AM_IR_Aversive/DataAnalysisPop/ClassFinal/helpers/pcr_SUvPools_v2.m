close all

whichClass   = 'Nspk';

crAmt = 0.001;

% Data settings
fn = set_paths_directories('','',1);

savedir = fullfile(fn.figs,'ClassResults',whichClass);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
% % 
% % datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass);
% % q = load(fullfile(fn.processed,'Units'));
% % UnitInfo_AM = q.UnitInfo;
% % clear q
% % 
% % % Load SU data 
% % q=load(fullfile(datadir_AM,'each','CR_each.mat'));
% % CReach_AM = q.CR;
% % clear q
% % 
% % % Load subpop results    
% % % q=load(fullfile(datadir_AM,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
% % % CR_pfr_AM = q.CR;
% % % clear q
% % q=load(fullfile(datadir_AM,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
% % CR_Qpfr_AM = q.CR;
% % clear q
% % q=load(fullfile(datadir_AM,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
% % CR_dp_AM = q.CR;
% % clear q
% % 
% % % Also get CTTS
% % [CTTS_AM,theseCells_AM] = recallDataParams('AC','each',12);
% % 

%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitInfo_Sp = q.UnitInfo;
clear q

% Load SU data 
q=load(fullfile(datadir_Sp,'each','CR_each.mat'));
CReach_Sp = q.CR;
clear q

% Load subpop results    
% q=load(fullfile(datadir_Sp,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
% CR_pfr_Sp = q.CR;
% clear q
q=load(fullfile(datadir_Sp,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
CR_Qpfr_Sp = q.CR;
clear q
q=load(fullfile(datadir_Sp,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
CR_dp_Sp = q.CR;
clear q

% Also get CTTS
[CTTS_Sp,theseCells_Sp] = recallDataParams('Speech','each',12);


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
sqsmall   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];


%==========================================================================
%%                               FIG 4
%                    pooled data classifier results
%==========================================================================


% RS indices
iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);

% Significance check
UnSig_AM = bootstrap4significance(CReach_AM(iRS_AM,:));
UnSig_Sp = bootstrap4significance(CReach_Sp(iRS_Sp,:));

% Sort units by peakFR and dprime
% AM
[~,ipkFR_AM]   = rankPeakFR(CTTS_AM(iRS_AM,:,:,:));
[~,iSUdps_AM]  = sort(CReach_AM(iRS_AM,:).dprime,'descend');
% Speech
[~,ipkFR_Sp]   = rankPeakFR(CTTS_Sp(iRS_Sp,:,:,:));
[~,iSUdps_Sp]  = sort(CReach_Sp(iRS_Sp,:).dprime,'descend');


ymaxval = 5;

% Plot
hf=figure; 
set(hf,'Position',sqsmall)

%% ~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

% Peak FR  ::  Quantiles

% NC_pfr_AM = mode(CR_pfr_AM.nC);
% iCRs_pfr_AM = find(CR_pfr_AM.nC==NC_pfr_AM);
NC_dp_AM = mode(CR_dp_AM.nC);
iCRs_dp_AM = find(CR_dp_AM.nC==NC_dp_AM);

% iSig_AM = find(UnSig_AM(iSUdps_AM));
% iNS_AM  = find(~UnSig_AM(iSUdps_AM));
%         plot( find( ~UnSig(isort) ), dps( isort( ~UnSig(isort) ) ) )

subplot(2,2,3)
hold on
plot( find(~UnSig_AM(ipkFR_AM)),CReach_AM(iRS_AM,:).dprime(ipkFR_AM(~UnSig_AM(ipkFR_AM))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_AM(ipkFR_AM)),CReach_AM(iRS_AM,:).dprime(ipkFR_AM(UnSig_AM(ipkFR_AM))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% plot(1:length(ipkFR_AM),CReach_AM(iRS_AM(ipkFR_AM),:).dprime,'.k','MarkerSize',15)
% plot(CR_pfr_AM.iC(iCRs_pfr_AM),CR_pfr_AM.dprime(iCRs_pfr_AM),'.m','MarkerSize',25)

plot( (CR_Qpfr_AM.iC-[zeros(size(CR_Qpfr_AM.iC,1),1)  CR_Qpfr_AM.nC-1])', [CR_Qpfr_AM.dprime CR_Qpfr_AM.dprime]',...
    '-g','MarkerSize',25,'LineWidth',3)


xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(ipkFR_AM),6)))
xlim([0 numel(iRS_AM)])
ylim([-0.1 ymaxval])
title(['AM, pkFR'])%, NC=' num2str(NC_pfr_AM)])


% dprime  ::  vary pool size

subplot(2,2,1)
hold on
plot( find(~UnSig_AM(iSUdps_AM)),CReach_AM(iRS_AM,:).dprime(iSUdps_AM(~UnSig_AM(iSUdps_AM))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_AM(iSUdps_AM)),CReach_AM(iRS_AM,:).dprime(iSUdps_AM(UnSig_AM(iSUdps_AM))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% plot(1:length(iSUdps_AM),CReach_AM(iRS_AM(iSUdps_AM),:).dprime,'.k','MarkerSize',15)
% plot(CR_dp_AM.iC(iCRs_dp_AM),CR_dp_AM.dprime(iCRs_dp_AM),'.m','MarkerSize',25)

xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps_AM),6)))
xlim([0 numel(iRS_AM)])
ylim([-0.1 ymaxval])
title(['AM, best d'', NC=' num2str(NC_dp_AM)])


%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

% Peak FR  ::  Quantiles

% NC_pfr_Sp = mode(CR_pfr_Sp.nC);
% iCRs_pfr_Sp = find(CR_pfr_Sp.nC==NC_pfr_Sp);
NC_dp_Sp = mode(CR_dp_Sp.nC);
iCRs_dp_Sp = find(CR_dp_Sp.nC==NC_dp_Sp);

subplot(2,2,4)
hold on
plot( find(~UnSig_Sp(ipkFR_Sp)),CReach_Sp(iRS_Sp,:).dprime(ipkFR_Sp(~UnSig_Sp(ipkFR_Sp))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_Sp(ipkFR_Sp)),CReach_Sp(iRS_Sp,:).dprime(ipkFR_Sp(UnSig_Sp(ipkFR_Sp))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% plot(1:length(ipkFR_Sp),CReach_Sp(iRS_Sp(ipkFR_Sp),:).dprime,'.k','MarkerSize',15)
% plot(CR_pfr_Sp.iC(iCRs_pfr_Sp),CR_pfr_Sp.dprime(iCRs_pfr_Sp),'.m','MarkerSize',25)
plot( (CR_Qpfr_Sp.iC-[zeros(size(CR_Qpfr_Sp.iC,1),1)  CR_Qpfr_Sp.nC-1])', [CR_Qpfr_Sp.dprime CR_Qpfr_Sp.dprime]',...
    '-g','MarkerSize',25,'LineWidth',3)

xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(ipkFR_Sp),6)))
xlim([0 numel(iRS_Sp)])
ylim([-0.1 ymaxval])
title(['Speech, pkFR'])%, NC=' num2str(NC_pfr_Sp)])


% dprime  ::  vary pool size

subplot(2,2,2)
hold on
plot( find(~UnSig_Sp(iSUdps_Sp)),CReach_Sp(iRS_Sp,:).dprime(iSUdps_Sp(~UnSig_Sp(iSUdps_Sp))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_Sp(iSUdps_Sp)),CReach_Sp(iRS_Sp,:).dprime(iSUdps_Sp(UnSig_Sp(iSUdps_Sp))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% plot(1:length(iSUdps_Sp),CReach_Sp(iRS_Sp(iSUdps_Sp),:).dprime,'.k','MarkerSize',15)
% plot(CR_dp_Sp.iC(iCRs_dp_Sp),CR_dp_Sp.dprime(iCRs_dp_Sp),'.m','MarkerSize',25)

xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps_Sp),6)))
xlim([0 numel(iRS_Sp)])
ylim([-0.1 ymaxval])
title(['Speech, best d'', NC=' num2str(NC_dp_Sp)])


%========================
% Add varied pool size
%========================

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

N = histcounts(CR_dp_AM.iC,'BinMethod','integers');
iCs_dp_AM = find(N>4);

subplot(2,2,1)
hold on
for ic = 1:numel(iCs_dp_AM)
    CRidx = find(ismember(CR_dp_AM.iC,iCs_dp_AM(ic)));
    [ncs,incs] = sort(CR_dp_AM(CRidx,:).nC);
    
    plot(iCs_dp_AM(ic)-1+ncs,CR_dp_AM.dprime(CRidx(incs)),'-b','LineWidth',3)
end


%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

N = histcounts(CR_dp_Sp.iC,'BinMethod','integers');
iCs_dp_Sp = find(N>4);

subplot(2,2,2)
hold on
for ic = 1:numel(iCs_dp_Sp)
    CRidx = find(ismember(CR_dp_Sp.iC,iCs_dp_Sp(ic)));
    [ncs,incs] = sort(CR_dp_Sp(CRidx,:).nC);
    
    plot(iCs_dp_Sp(ic)-1+ncs,CR_dp_Sp.dprime(CRidx(incs)),'-b','LineWidth',3)
end


% keyboard
print_eps_kp(gcf,fullfile(savedir,'SUvPools_all_5'))




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


%%


