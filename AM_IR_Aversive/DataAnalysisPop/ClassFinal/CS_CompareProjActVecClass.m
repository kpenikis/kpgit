% compare ActVec and Projections classifiers

% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

%%
% results for adding RS cells dp ranked

close all

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/allRS/CR_vFull_allRS.mat')
CR_all_Proj = CR;
clear CR
CR_all_Proj = CR_all_Proj(1,:);

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/dpRank_RS/CR_vFull_dpRank_RS.mat')
CR_Proj = CR;
clear CR

CR_Proj = CR_Proj(CR_Proj.iC==1,:);
[nC_P,inC_P] = sort(CR_Proj.nC);
CR_Proj = CR_Proj(inC_P,:);


load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/ActVec/allRS/CR_vActVec_allRS.mat')
CR_all_ActV = CR;
clear CR

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/ActVec/dpRank_RS/CR_vActVec_dpRank_RS.mat')
CR_ActV = CR;
clear CR

CR_ActV = CR_ActV(CR_ActV.iC==1,:);
[nC_A,inC_A] = sort(CR_ActV.nC);
CR_ActV = CR_ActV(inC_A,:);



hf=figure;

subplot(2,1,1)
ip(1)=plot([nC_P; 181],[CR_Proj.PC; CR_all_Proj.PC],'k','LineWidth',3);
hold on
ip(2)=plot([nC_A; 181],[CR_ActV.PC; mean(CR_all_ActV.PC)],'b','LineWidth',3);
legend(ip,{'Projections' 'ActVec'},'Location','best')
% xlabel('# Cells in ensemble, added from best to worst SU d''')
ylabel('Percent Correct')
set(gca,'Color','none')
xlim([0 181])
ylim([60 100])


subplot(2,1,2)
ip(1)=plot([nC_P; 181],[CR_Proj.dprime; CR_all_Proj.dprime],'k','LineWidth',3);
hold on
ip(2)=plot([nC_A; 181],[CR_ActV.dprime; mean(CR_all_ActV.dprime)],'b','LineWidth',3);
legend(ip,{'Projections' 'ActVec'},'Location','best')
xlabel('# Cells in ensemble, added from best to worst SU d''')
ylabel('d''')
set(gca,'Color','none')
xlim([0 181])
ylim([0 5])


fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults','Corrections');
print_eps_kp(hf,fullfile(savedir,'ProjActVec_addCells'))



%%
% results for quantiles

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/Q_pkFR/CR_vFull_Q_pkFR.mat')
CR_Proj = CR;
clear CR

[iC_P,iiC_P] = sort(CR_Proj.iC);
CR_Proj = CR_Proj(iiC_P,:);


load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/ActVec/Q_pkFR/CR_vActVec_Q_pkFR.mat')
CR_ActV = CR;
clear CR

[iC_A,iiC_A] = sort(CR_ActV.iC);
CR_ActV = CR_ActV(iiC_A,:);


subplot(2,2,3)

ip(1)=plot(iC_P,CR_Proj.dprime,'k','LineWidth',2);
hold on
ip(2)=plot(iC_A,CR_ActV.dprime,'b','LineWidth',2);
legend(ip,{'Projections' 'ActVec'},'Location','best')
xlabel('Quantiles, 36 cells starting with N')
ylabel('d''')


%%
% results for random 10 cells

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/Rand_RS/CR_Full_Rand_RS.mat')
CR_Proj = CR;
clear CR

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/ActVec/Rand_RS/CR_ActVec_Rand_RS.mat')
CR_ActV = CR;
clear CR

subplot(2,2,4)
plotSpread(CR_Proj.dprime,'distributionIdx',ones(size(CR_Proj.dprime)),'distributionColors','k','showMM',1)
hold on
plotSpread(CR_ActV.dprime,'distributionIdx',1+ones(size(CR_AV.dprime)),'distributionColors','b','showMM',1)
