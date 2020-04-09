close all

whichClass   = 'Full';

crAmt = 0.001;
nStim = 8;

% Data settings
fn = set_paths_directories('','',1);

savedir = fullfile(fn.figs,'ClassResults',whichClass);
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
% q=load(fullfile(datadir_AM,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
% CR_pfr_AM = q.CR;
% clear q
% q=load(fullfile(datadir_AM,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
% CR_Qpfr_AM = q.CR;
% clear q
% q=load(fullfile(datadir_AM,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
% CR_dp_AM = q.CR;
% clear q

% Also get CTTS
[CTTS_AM,theseCells_AM] = recallDataParams('AC','each',12);


%~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~

% datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
% q = load(fullfile(fn.processed,'UnitsVS'));
% UnitInfo_Sp = q.UnitInfo;
% UnitData_Sp = q.UnitData;
% clear q
% 
% % Load SU data 
% q=load(fullfile(datadir_Sp,'each','CR_each.mat'));
% CReach_Sp = q.CR;
% clear q
% 
% % Load subpop results    
% q=load(fullfile(datadir_Sp,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
% CR_pfr_Sp = q.CR;
% clear q
% q=load(fullfile(datadir_Sp,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
% CR_Qpfr_Sp = q.CR;
% clear q
% q=load(fullfile(datadir_Sp,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
% CR_dp_Sp = q.CR;
% clear q
% 
% % Also get CTTS
% [CTTS_Sp,theseCells_Sp] = recallDataParams('Speech','each',12);


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall   = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
sqsmall     = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];


%% ~~~~ d'/PC vals for each stimulus

% CellTypes
iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
iNS = find(UnitInfo_AM(theseCells_AM,:).TroughPeak<=0.43);


[dps_AM,iSUdps_AM] = sort(CReach_AM(iRS_AM,:).dprime,'descend');

% 
% iRS(iSUdps_AM)
% 
% maxdps = max(dpStim,[],2);
% [maxdps_pl,imaxdps] = sort(maxdps,'descend');
% 
% figure;
% plot(1:length(dps_AM),dps_AM,'.k','MarkerSize',15)
% hold on
% plot(1:numel(maxdps_pl),maxdps_pl,'.m','MarkerSize',15)
% xlim([0 length(dps_AM)])
% ylim([-0.1 5])
% set(gca,'Color','none')
% grid on
% 
% print_eps_kp(gcf,fullfile(figsavedir,'SUdps_sort'))



% figure;
% % set(gcf,'Position',widesmall)
% plot(1:length(dps),dps,'.k','MarkerSize',15)
% xlim([0 length(dps)])
% ylim([-0.05 2.5])
% set(gca,'ytick',0:0.5:2.5,'xtick',[5 10 30 50 100 150])
% grid on
% title([whichStim ' -- SU dprimes'])
% ylabel('d prime')
% xlabel('Cell number')

% keyboard
% print_eps_kp(gcf,fullfile(savedir,'cdf_SUdps'))


% For each cell, get its PC & d' for each stimulus
pcStim_AM = nan(length(iRS_AM),nStim);
dpStim_AM = nan(length(iRS_AM),nStim);
Sparss_AM = nan(length(iRS_AM),1);
for ii = 1:length(iRS_AM)
        
    ConfMat = mean(CReach_AM(iRS_AM(ii),:).Results{:},3,'omitnan');
    
    pcStim_AM(ii,:) = diag(ConfMat)';
    
    ConfMat(ConfMat==0) = 0.01;
    ConfMat(ConfMat==1) = 0.99;
    for ist = 1:size(ConfMat,1)
        othSt = 1:size(ConfMat,1);
        dpStim_AM(ii,ist) =  norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(othSt~=ist,ist)),0,1);
    end
    Sparss_AM(ii) = calculateSparseness(dpStim_AM(ii,:)');
end
% Cell rank is a little different when d' calculated this way
% [meandps,iUdp] = sort(mean(dpStim_AM,2),'descend');

UnSig_AM = bootstrap4significance(CReach_AM(iRS_AM,:));

avg_dp = CReach_AM(iRS_AM,:).dprime;



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%% Plot mean dprime vs sparseness of performance across stimuli
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure;
set(gcf,'Position',sqsmall)
subplot(1,4,1:3)
plot(avg_dp(~UnSig_AM),Sparss_AM(~UnSig_AM),'xk')
hold on
plot(avg_dp(UnSig_AM),Sparss_AM(UnSig_AM),'.','MarkerSize',20)
xlabel('mean d'', each SU')
ylabel('stim sparseness')
ylim([0 1])
xlim([-0.2 2.5])
axis square

subplot(1,4,4)
% histogram(Sparss_AM(~UnSig_AM),0:0.05:1,'FaceColor','k','EdgeColor','k')
% hold on
histogram(Sparss_AM(UnSig_AM),0:0.05:1)
xlim([0 1])
view(90,-90)

print_eps_kp(gcf,fullfile(savedir,'SU_meanDP-spsns_AM'))



%% PLAY WITH 2 EXAMPLE CELLS
% same mean d' different sparseness

exC_data = find(abs(avg_dp-0.36673)<0.0001); %index out of 181
exC_iUn  = theseCells_AM(iRS_AM((abs(avg_dp-0.36673)<0.0001))); %index out of 277

Sparss_AM(exC_data)

UnitData_AM(exC_iUn(1))
% mean d' = 0.37
% st sps  = 0.24

UnitData_AM(exC_iUn(2))
% mean d' = 0.37
% st sps  = 0.73

figure;
plot(dpStim_AM(exC_data,:)','LineWidth',3)
xlabel('stim')
ylabel('d''')
print_eps_kp(gcf,fullfile(savedir,'SU_exCs_samemean_stimtuning_AM'))

% Load rand subpop results    
q=load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/PoolRand/RS/CR_vPoolRand_RS.mat');
CR_rand_AM = q.CR;
clear q

idx = iRS_AM(exC_data); %index out of 235
dp_c1 = [];
dp_c2 = [];
for ii = 1:size(CR_rand_AM,1)
    if ismember(idx(1),CR_rand_AM.SUids{ii})
        dp_c1 = [dp_c1; CR_rand_AM.dprime(ii)];
    end
    if ismember(idx(2),CR_rand_AM.SUids{ii})
        dp_c2 = [dp_c2; CR_rand_AM.dprime(ii)];
    end
end

figure;
plotSpread([dp_c1; dp_c2],'distributionIdx',[ones(size(dp_c1)); 1+ones(size(dp_c2))],'distributionColors','k')
hold on
boxplot([dp_c1; dp_c2],[ones(size(dp_c1)); 1+ones(size(dp_c2))])
ylim([-0.2 2.5])
ylabel('d''')
xlabel('ex cells, same mean d''')
set(gca,'xticklabel',{'low sps' 'high sps'})
print_eps_kp(gcf,fullfile(savedir,'SU_exCs_samemean_randpools_AM'))



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%% Sort by mean dprime, plot min to max across stimuli
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

figure;
set(gcf,'Position',sqsmall)
plot([avg_dp() avg_dp()]'...
    ,[min(dpStim_AM(),[],2) max(dpStim_AM(),[],2)]','-k')
hold on
plot(avg_dp(), avg_dp,'.k','MarkerSize',10)
plot(avg_dp(), min(dpStim_AM(),[],2),'.m','MarkerSize',20)
plot(avg_dp(), max(dpStim_AM(),[],2),'.g','MarkerSize',20)
xlabel('mean d'', each SU')
ylabel('max (pink) and min (green) d'' across stimuli')
xlim([-0.2 2.5])
axis square

print_eps_kp(gcf,fullfile(savedir,'SU_meanDP-maxmin_AM'))

% Same, in log space
figure;
set(gcf,'Position',sqsmall)
plot([avg_dp() avg_dp()]'...
    ,[min(dpStim_AM(),[],2) max(dpStim_AM(),[],2)]','-k')
hold on
plot(avg_dp(), avg_dp,'.k','MarkerSize',10)
plot(avg_dp(), min(dpStim_AM(),[],2),'.m','MarkerSize',20)
plot(avg_dp(), max(dpStim_AM(),[],2),'.g','MarkerSize',20)
xlabel('mean d'', each SU')
ylabel('max (pink) and min (green) d'' across stimuli')
set(gca,'xscale','log','yscale','log')
xlim(10.^[-3 1])
ylim(10.^[-3 1])
axis square

print_eps_kp(gcf,fullfile(savedir,'SU_log_meanDP-maxmin_AM'))



%% Sort by mean d', plot d' of each stimulus in subplots

StimOrder_AM = [8 3 2 5 6 7 4 1];
StimOrder_Sp = [4 3 6 1 2 5 7 8];
newcolors = cmocean('thermal',8);


hfa=figure;
set(gcf,'Position',fullscreen./[1 1 1 2])
    
for ist = StimOrder_AM
    
%     hf(ist)=figure;
%     set(gcf,'Position',widesmall./[1 1 2 1])
    subplot(2,4,ist)
    hold on
    
    % Add this stimulus dp for each cell
    plot([1:length(dps_AM); 1:length(dps_AM)],[zeros(length(iSUdps_AM),1) dpStim_AM(iSUdps_AM,ist)]',...
        '-','Color',newcolors(ist==StimOrder_AM,:),'LineWidth',2)
    
    plot(1:length(dps_AM),dps_AM,'.k','MarkerSize',10)
    
    xlim([0 length(dps_AM)+1])
    ylim([-0.05 2.5])
    set(gca,'ytick',0:0.5:2.5,'xtick',[0 length(iSUdps_AM)])%,'Color','none')
    grid on
    title(['SU dprimes -- st# ' num2str(ist)])
    if ist==1
        ylabel('d prime')
    end    
    if ist==5
        ylabel('d prime')
        xlabel('Cell number')
    end
    
%     print_eps_kp(hf(ist),fullfile(savedir,['SUdps_AM_st' num2str(ist)]))
end

print_eps_kp(hfa,fullfile(savedir,'SUdps_AM_all'))




keyboard


%% Compare d' to FF 







%%
% nSpks  = (mean(sum(mean(CTTS,3,'omitnan'),2,'omitnan'),4));
% nSpkRS = (mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2,'omitnan'),4));
% nSpkNS = (mean(sum(mean(CTTS(iNS,:,:,:),3,'omitnan'),2,'omitnan'),4));


figure;
subplot(1,4,1)
hold on

% Manually make boxplots
q5  = quantile(CReach.dprime,0.05);
q25 = quantile(CReach.dprime,0.25);
q75 = quantile(CReach.dprime,0.75);
q95 = quantile(CReach.dprime,0.95);

plot([1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
fill(1+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')

ylabel('dprime')
set(gca,'Color','none')
box off
ylim([-0.5 2.5])
plotSpread(CReach.dprime,'distributionIdx',ones(size(CReach.dprime)),'distributionColors','k','showMM',3)


subplot(1,4,2:4)
%         plot(nSpks,CReach.dprime,'k.')
hold on
plot(nSpkRS,CReach.dprime(iRS),'r.')
plot(nSpkNS,CReach.dprime(iNS),'b.')
ylim([-0.5 2.5])
set(gca,'ytick',[],'Color','none')
box off
xlabel('N spks')
%         title('SU dprimes vs mean N spks per stim')


% Stats
[r_RS,p_RS] = corr(nSpkRS,CReach.dprime(iRS),'Type','Spearman');
[r_NS,p_NS] = corr(nSpkNS,CReach.dprime(iNS),'Type','Spearman');
[r,p]       = corr(nSpks,CReach.dprime,'Type','Spearman');

text(40,-0.3,sprintf('Spearman r=%0.2f, p=%0.3e',r,p))


savename = sprintf('SU_dps_%s_%s',varPar,whichStim);

keyboard

print_eps_kp(gcf,fullfile(savedir,savename))

