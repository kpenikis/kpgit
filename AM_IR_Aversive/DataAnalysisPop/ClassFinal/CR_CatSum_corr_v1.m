function CR_CatSum_corr_v1
% Plots Fig 10* last panel
% new data

% close all

whichClass   = 'Full';
whichCells   = 'dpRank_RS';

whichStim    = 'Speech';

minNtrs      = 17;

dpmaxval     = 4;
crAmt        = 0.001;

% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults','Corrections','CatSum');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


switch whichStim
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'AC'
        
        rootdir = fullfile(fn.figs,'ClassAM',whichStim,whichClass,['minTrs' num2str(minNtrs)]);
        switch whichClass
            case 'Full'
                rootdir_sum = fullfile(fn.figs,'ClassAM',whichStim,'Sum',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'Sum';
            case 'ActVec'
                rootdir_sum = fullfile(fn.figs,'ClassAM',whichStim,'SumAV',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'SumAV';
        end
        
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'Speech'
        
        rootdir = fullfile(fn.figs,'ClassSpeech',whichStim,whichClass,['minTrs' num2str(minNtrs)]);
        switch whichClass
            case 'Full'
                rootdir_sum = fullfile(fn.figs,'ClassSpeech',whichStim,'Sum',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'Sum';
            case 'ActVec'
                rootdir_sum = fullfile(fn.figs,'ClassSpeech',whichStim,'SumAV',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'SumAV';
        end
        
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
end


% Load SU data 
q=load(fullfile(rootdir,'CR_each.mat'));
CReach = q.CR;
clear q

% Load subpop results    
% Cat
q=load(fullfile(rootdir,whichCells,['CR_v' whichClass '_' whichCells '.mat']));
CR_Cat = q.CR;
clear q
% Sum
q=load(fullfile(rootdir_sum,whichCells,['CR_v' sumclassstr  '_' whichCells '.mat']));
CR_Sum = q.CR;
clear q


% Also get CTTS
[~,theseCells] = recallDataParams(whichStim,'each',minNtrs);


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen  = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
sqsmall   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];


%==========================================================================
%%                               FIG X
%                    pooled data classifier results
%==========================================================================

nStim    = size(CReach(1,:).Results{:},1);

% RS indices
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);

% Significance check
% UnSig_AM = bootstrap4significance(CReach_AM(iRS_AM,:));
% UnSig_Sp = bootstrap4significance(CReach_Sp(iRS_Sp,:));

% Sort units by peakFR and dprime
[~,iSUdps]  = sort(CReach(iRS,:).dprime,'descend');

% Filter Cat data to just FC = 1
CR_Cat = CR_Cat(ismember(CR_Cat.iC,1),:);


%% Plot

hf=figure; 
set(hf,'Position',sqsmall)


% ~~~~~~~~~~~~~~~~~~~~~~~~~ SU data

subplot(2,2,1)
hold on
% plot( find(~UnSig_AM(iSUdps_AM)),CReach_AM(iRS_AM,:).dprime(iSUdps_AM(~UnSig_AM(iSUdps_AM))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
% plot( find(UnSig_AM(iSUdps_AM)),CReach_AM(iRS_AM,:).dprime(iSUdps_AM(UnSig_AM(iSUdps_AM))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
dp = CReach(iRS(iSUdps),:).dprime;
dp(dp>dpmaxval) = dpmaxval;
plot(1:length(iSUdps),dp,'.k','MarkerSize',15)


% ~~~~~~~~~~~~~~~~~~~~~~~~~ Pooled 

% Cat
[ncs,incs] = sort(CR_Cat.nC);
dp = CR_Cat.dprime(incs);
dp(dp>dpmaxval) = dpmaxval;
plot(ncs,dp,':r','LineWidth',3)

% Sum
[ncs,incs] = sort(CR_Sum.nC);
dp = CR_Sum.dprime(incs);
dp(dp>dpmaxval) = dpmaxval;
plot(ncs,dp,'-r','LineWidth',3)

% ~~~~~~~~~~~~~~~~~~~~~~~~~ Finish plot
xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps),6)))
xlim([0 numel(iRS)])
ylim([-0.3 dpmaxval])
title([whichStim ', rank best d'''])



%% Collect d' for each stimulus

% - - - RS - - - 
% SU
dp_RS_SU = nan(size(iSUdps,1),nStim);
for ii = 1:numel(iSUdps)
    dp_RS_SU(ii,:) = dp_from_ConfMat(mean(CReach(iRS(iSUdps(ii)),:).Results{:},3,'omitnan'),crAmt);
end

% Cat
dp_RS_Cat  = nan(size(CR_Cat,1),nStim);
Spsns_Cat  = nan(size(CR_Cat,1),1);
NC_Cat     = nan(size(CR_Cat,1),1);
[ncs,incs] = sort(CR_Cat.nC);
for ii = 1:numel(incs)
    dp_RS_Cat(ii,:) = dp_from_ConfMat(mean(CR_Cat(incs(ii),:).Results{:},3,'omitnan'),crAmt);
    Spsns_Cat(ii)   = calculateSparseness(dp_RS_Cat(ii,:));
    NC_Cat(ii)      = ncs(ii);
end
% dp_RS_Cat(end,:)    = dp_from_ConfMat(mean(CR_Cat.Results{:},3,'omitnan'),crAmt);
% Spsns_Cat(end)      = calculateSparseness(dp_RS_Cat(end,:));
% NC_Cat(end)         = numel(CR_Cat.SUdps{:});

% Sum
dp_RS_Sum  = nan(size(CR_Sum,1),nStim);
Spsns_Sum  = nan(size(CR_Sum,1),1);
NC_Sum     = nan(size(CR_Sum,1),1);
[ncs,incs] = sort(CR_Sum.nC);
for ii = 1:numel(incs)
    dp_RS_Sum(ii,:) = dp_from_ConfMat(mean(CR_Sum(incs(ii),:).Results{:},3,'omitnan'),crAmt);
    Spsns_Sum(ii)   = calculateSparseness(dp_RS_Sum(ii,:));
    NC_Sum(ii)      = ncs(ii);
end
% dp_RS_Sum(end,:)    = dp_from_ConfMat(mean(CR_Sum.Results{:},3,'omitnan'),crAmt);
% Spsns_Sum(end)      = calculateSparseness(dp_RS_Sum(end,:));
% NC_Sum(end)         = numel(CR_Sum.SUdps{:});


% Correct each to max d' value
dp_RS_SU(dp_RS_SU>dpmaxval)   = dpmaxval;
dp_RS_Cat(dp_RS_Cat>dpmaxval) = dpmaxval;
dp_RS_Sum(dp_RS_Sum>dpmaxval) = dpmaxval;


%% Add sparseness to task avg 

% % % Cat
% % plot(NC_Cat,Spsns_Cat,':k','LineWidth',3)
% % plot(NC_Sum,Spsns_Sum,'-k','LineWidth',3)
% % 
% % plot(NC_Cat,min(dp_RS_Cat,[],2),':b','LineWidth',3)
% % plot(NC_Sum,min(dp_RS_Sum,[],2),'-b','LineWidth',3)
% % 
% % 
% % 
% % % Cat
% % [bestdp_Cat,imax] = max(CR_dp_Cat.dprime);
% % bestNC_Cat = CR_dp_Cat.nC(imax);
% % mindp_Cat  = min(dp_RS_Cat(bestNC_Cat==NC_Cat,:));
% % 
% % % Sum 
% % [bestdp_Sum,imax] = max(CR_dp_Sum.dprime);
% % bestNC_Sum = CR_dp_Sum.nC(imax);
% % mindp_Sum  = min(dp_RS_Sum(bestNC_Sum==NC_Sum,:));
% % 
% % figure;
% % plot([0 5],[0 5],'-k')
% % hold on
% % plot(mindp_Cat,bestdp_Cat,'og')
% % plot(min(dp_RS_Cat(end,:)),CR_Cat.dprime,'xg')
% % % bestdp_Cat - mindp_Cat
% % % CR_Cat.dprime - min(dp_RS_Cat(end,:))
% % plot(mindp_Sum,bestdp_Sum,'om')
% % plot(min(dp_RS_Sum(end,:)),CR_Sum.dprime,'xm')
% % axis square


%% Add d' with all cells, each stimulus

subplot(2,2,3)
hold on

% set stim order
diffs = dp_RS_Sum(end,:) - dp_RS_Cat(end,:);
[dpd,idpd] = sort(diffs,'descend');

plot(1:nStim,dp_RS_Cat(end,idpd),':r','LineWidth',3)
plot(1:nStim,dp_RS_Sum(end,idpd),'-r','LineWidth',3)

% ~~~~~~~~~~~~~~~~~~~~~~~~~ Finish plot
xlabel('Stimulus')
ylabel('d'' all RS cells')
grid off
box off
set(gca,'Color','none','xtick',1:nStim,'xticklabel',idpd)
ylim([0 dpmaxval])
xlim([0.5 nStim+0.5])


subplot(2,2,4)
bar(1:nStim,dpd,'r')

% ~~~~~~~~~~~~~~~~~~~~~~~~~ Finish plot
xlabel('Stimulus')
ylabel('d'' Sum - d'' Cat')
grid off
box off
set(gca,'Color','none','xtick',1:nStim,'xticklabel',idpd)
ylim([-3.25 1.333])
% ylim([-1.5 1])
xlim([0.5 nStim+0.5])


savename = ['AddCells_' whichStim '_' whichClass '_' num2str(minNtrs) 'trs'];
print_eps_kp(hf,fullfile(savedir,savename))

return




%% Plot add cells, each stimulus

hf_supp=figure; 
set(hf_supp,'Position',halfscreen)

for ist = 1:nStim
    
    subplot(4,2,ist)
    hold on
    
    plot(1:size(dp_RS_SU,1),dp_RS_SU(:,ist),'.k','MarkerSize',10)
    
    % Cat
    plot(NC_Cat,dp_RS_Cat(:,ist),':','Color','r','LineWidth',3)
    
    % Sum
    plot(NC_Sum,dp_RS_Sum(:,ist),'-r','LineWidth',3)
    
    
%     xlabel('Cell N')
    ylabel('d''')
    grid on
    box off
    set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps),6)))
    xlim([0 numel(iRS)])
    ylim([-0.3 dpmaxval+0.5])
    title(ist)
end

print_eps_kp(hf_supp,fullfile(savedir,[savename '_byStim']))


%%

return








%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

% Peak FR  ::  Quantiles

% NC_pfr_Sp = mode(CR_pfr_Sp.nC);
% iCRs_pfr_Sp = find(CR_pfr_Sp.nC==NC_pfr_Sp);
% NC_dp_Sp = mode(CR_dp_Sp.nC);
% iCRs_dp_Sp = find(CR_dp_Sp.nC==NC_dp_Sp);
% 
% subplot(2,2,4)
% hold on
% plot( find(~UnSig_Sp(ipkFR_Sp)),CReach_Sp(iRS_Sp,:).dprime(ipkFR_Sp(~UnSig_Sp(ipkFR_Sp))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
% plot( find(UnSig_Sp(ipkFR_Sp)),CReach_Sp(iRS_Sp,:).dprime(ipkFR_Sp(UnSig_Sp(ipkFR_Sp))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% % plot(1:length(ipkFR_Sp),CReach_Sp(iRS_Sp(ipkFR_Sp),:).dprime,'.k','MarkerSize',15)
% % plot(CR_pfr_Sp.iC(iCRs_pfr_Sp),CR_pfr_Sp.dprime(iCRs_pfr_Sp),'.m','MarkerSize',25)
% plot( (CR_Qpfr_Sp.iC-[zeros(size(CR_Qpfr_Sp.iC,1),1)  CR_Qpfr_Sp.nC-1])', [CR_Qpfr_Sp.dprime CR_Qpfr_Sp.dprime]',...
%     '-g','MarkerSize',25,'LineWidth',3)
% 
% xlabel('Cell N')
% ylabel('d'' across stimuli')
% grid on
% box off
% set(gca,'Color','none','xtick',round(linspace(0,length(ipkFR_Sp),6)))
% xlim([0 numel(iRS_Sp)])
% ylim([-0.1 ymaxval])
% title(['Speech, pkFR'])%, NC=' num2str(NC_pfr_Sp)])
% 
% 
% % dprime  ::  SU data
% 
% subplot(2,2,2)
% hold on
% plot( find(~UnSig_Sp(iSUdps_Sp)),CReach_Sp(iRS_Sp,:).dprime(iSUdps_Sp(~UnSig_Sp(iSUdps_Sp))),'+','Color',0.5*[1 1 1],'MarkerSize',5)
% plot( find(UnSig_Sp(iSUdps_Sp)),CReach_Sp(iRS_Sp,:).dprime(iSUdps_Sp(UnSig_Sp(iSUdps_Sp))),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])
% % plot(1:length(iSUdps_Sp),CReach_Sp(iRS_Sp(iSUdps_Sp),:).dprime,'.k','MarkerSize',15)
% % plot(CR_dp_Sp.iC(iCRs_dp_Sp),CR_dp_Sp.dprime(iCRs_dp_Sp),'.m','MarkerSize',25)
% 
% xlabel('Cell N')
% ylabel('d'' across stimuli')
% grid on
% box off
% set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps_Sp),6)))
% xlim([0 numel(iRS_Sp)])
% ylim([-0.1 ymaxval])
% title(['Speech, best d'', NC=' num2str(NC_dp_Sp)])


%=========================
% Add pool data (rank d')
%=========================

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

% N = histcounts(CR_dp_Sp.iC,'BinMethod','integers');
% iCs_dp_Sp = find(N>4);
% 
% subplot(2,2,2)
% hold on
% for ic = 1:numel(iCs_dp_Sp)
%     CRidx = find(ismember(CR_dp_Sp.iC,iCs_dp_Sp(ic)));
%     [ncs,incs] = sort(CR_dp_Sp(CRidx,:).nC);
%     
%     plot(iCs_dp_Sp(ic)-1+ncs,CR_dp_Sp.dprime(CRidx(incs)),'-b','LineWidth',3)
% end


% keyboard
% print_eps_kp(gcf,fullfile(savedir,'SUvPools_all_5'))



%==========================================================================
%% FIG 3 
%  sort d', cumulative N spikes
%==========================================================================


% Also get CTTS
[CTTS,theseCells] = recallDataParams(whichStim,'each',12);

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

[dps_AM,iSUdps] = sort(CReach(iRS,:).dprime,'ascend');


idp1 = find(dps_AM<prctile(dps_AM,94),1,'last');
avgNspks = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
pBest_nSpk = (sum(avgNspks(iSUdps((idp1+1):end)))/sum(avgNspks) + sum(avgNspks(iSUdps((idp1):end)))/sum(avgNspks)) /2;
pRest_nSpk = 1-pBest_nSpk;

dpmaxval = 2.5;
figure; 

yyaxis right
hold on
plot( find(~UnSig_AM(iSUdps)),dps_AM(~UnSig_AM(iSUdps)),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_AM(iSUdps)),dps_AM(UnSig_AM(iSUdps)),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])

plot([idp1 idp1]+0.5,[-0.1 dpmaxval],'-','Color',0.7*[1 1 1])
ylabel('dprime')
xlim([0 numel(dps_AM)+1])
ylim([-0.1 dpmaxval])

yyaxis left
plot(1:numel(dps_AM),cumsum(avgNspks(iSUdps))./sum(avgNspks(iSUdps)),'-','LineWidth',3)
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


dpmaxval = 3;
figure; 

yyaxis right
hold on
plot( find(~UnSig_Sp(iSUdps_Sp)),dps_Sp(~UnSig_Sp(iSUdps_Sp)),'+','Color',0.5*[1 1 1],'MarkerSize',5)
plot( find(UnSig_Sp(iSUdps_Sp)),dps_Sp(UnSig_Sp(iSUdps_Sp)),'o','MarkerEdgeColor',0.02*[1 1 1],'MarkerSize',5.5,'MarkerFaceColor',0.51*[1 1 1])

plot([idp1 idp1]+0.5,[-0.1 dpmaxval],'-','Color',0.7*[1 1 1])
ylabel('dprime')
xlim([0 numel(dps_Sp)+1])
ylim([-0.1 dpmaxval])

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
plot(CRe_Proj_AM.dprime,CReach.dprime,'.k')
hold on
plot([-0.2 3],[-0.2 3],'-')
axis square
set(gca,'Color','none','ytick',0:3,'xtick',0:3)
xlim([-0.2 3])
ylim([-0.2 3])

[r,p] = corr(CRe_Proj_AM.dprime,CReach.dprime,'type','Spearman')
p_wsr_AM = signrank(CRe_Proj_AM.dprime,CReach.dprime)
median(CReach.dprime-CRe_Proj_AM.dprime)

statstext = sprintf('SINUSOIDAL\nPearson r=%0.2f, p<0.001\nsignrank p=%0.2e\nd'' ActVec-Proj median = %0.2f',r,p_wsr_AM,median(CReach.dprime-CRe_Proj_AM.dprime));

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





end

