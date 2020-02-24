close all
crAmt = 0.001;

whichClass   = 'Full';
% whichCells   = 'pkFR_RS'; 'dpRank_RS';

% Data settings
fn = set_paths_directories('','',1);

% AM
datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass,'each');
q = load(fullfile(fn.processed,'Units'));
UnitInfo_AM = q.UnitInfo;
clear q

% Speech
datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass,'each');
q = load(fullfile(fn.processed,'UnitsVS'));
UnitInfo_Sp = q.UnitInfo;
clear q

% Load SU data
q=load(fullfile(datadir_AM,'CR_each.mat'));
CReach_AM = q.CR;
clear q
q=load(fullfile(datadir_Sp,'CR_each.mat'));
CReach_Sp = q.CR;
clear q


% Also get CTTS
[~,theseCells_AM] = recallDataParams('AC','each');
[~,theseCells_Sp] = recallDataParams('Speech','each');


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];


%% d' and sparseness, distributions for AM and Speech, iRS and iNS 

figure;
set(gcf,'Position',widesmall)

%%%%%%   AM

% Calculate basic reults measures
dpStim     = nan(8,size(CReach_AM,1));
Spsns_AM   = nan(1,size(CReach_AM,1));
Covrg_AM   = nan(1,size(CReach_AM,1));
for icr = 1:size(CReach_AM,1)
    ConfMat = mean(CReach_AM(icr,:).Results{:},3);
    dpStim(:,icr) = dp_from_ConfMat(ConfMat,0.001);
    Spsns_AM(icr) = calculateSparseness(dpStim(:,icr));
    Covrg_AM(icr) = 1-Spsns_AM(icr);
end

iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
iNS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak<=0.43);

subplot(1,2,1)
hold on
plotSpread(CReach_AM.dprime(iRS_AM),'distributionIdx',ones(1,length(iRS_AM)),'showMM',3)
plotSpread(CReach_AM.dprime(iNS_AM),'distributionIdx',4+ones(1,length(iNS_AM)),'showMM',3)
ylabel('d''')

subplot(1,2,2)
hold on
plotSpread(Covrg_AM(iRS_AM),'distributionIdx',ones(1,length(iRS_AM)),'showMM',3)
plotSpread(Covrg_AM(iNS_AM),'distributionIdx',4+ones(1,length(iNS_AM)),'showMM',3)
ylabel('Coverage')


%%%%%%   Speech

% Calculate basic reults measures
dpStim     = nan(8,size(CReach_Sp,1));
Spsns_Sp   = nan(1,size(CReach_Sp,1));
Covrg_Sp   = nan(1,size(CReach_Sp,1));
for icr = 1:size(CReach_Sp,1)
    ConfMat = mean(CReach_Sp(icr,:).Results{:},3);
    dpStim(:,icr) = dp_from_ConfMat(ConfMat,0.001);
    Spsns_Sp(icr) = calculateSparseness(dpStim(:,icr));
    Covrg_Sp(icr) = 1-Spsns_Sp(icr);
end

iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);
iNS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak<=0.43);

subplot(1,2,1)
hold on
plotSpread(CReach_Sp.dprime(iRS_Sp),'distributionColors','k','distributionIdx',1+ones(1,length(iRS_Sp)),'showMM',3)
plotSpread(CReach_Sp.dprime(iNS_Sp),'distributionColors','k','distributionIdx',5+ones(1,length(iNS_Sp)),'showMM',3)
ylabel('d''')
set(gca,'xtick',[1.5 5.5],'xticklabel',{'RS' 'NS'})

subplot(1,2,2)
hold on
plotSpread(Covrg_Sp(iRS_Sp),'distributionColors','k','distributionIdx',1+ones(1,length(iRS_Sp)),'showMM',3)
plotSpread(Covrg_Sp(iNS_Sp),'distributionColors','k','distributionIdx',5+ones(1,length(iNS_Sp)),'showMM',3)
ylabel('Coverage')
set(gca,'xtick',[1.5 5.5],'xticklabel',{'RS' 'NS'})

title(whichClass)

savedir = fullfile(fn.figs,'ClassResults',whichClass);
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(gcf,fullfile(savedir,'SU_AM-Sp_RSNS_comparison'))



keyboard


