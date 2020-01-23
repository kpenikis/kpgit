
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/PoolRand/RS/CR_vPoolRand_RS.mat')

% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');


figure;
plotSpread(CR.dprime,'showMM',1)
title('Random pools of 10 cells')
xlabel([num2str(round(sum(CR.dprime>1)/size(CR,1)*100)) ...
    '% of pools d''>1, '...
    num2str(round(sum(CR.dprime>2)/size(CR,1)*100)) ...
    '% d''>2, '])

print_eps_kp(gcf,fullfile('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/PoolRand/RS','RandPoolsDistr'))



figure;
% set(gcf,'Position',widesmall)
hold on
for icr = 1:size(CR,1)
    SUdps = CR.SUdps(icr);
    plotSpread(SUdps{:},'distributionIdx', CR.iC(icr)*ones(size(SUdps{:})),'showMM',3)
end

plot(CR.iC,CR.dprime,'k-')


figure; 
% set(gcf,'Position',widesmall)
hold on
for icr = 1:size(CR,1)
    SUdps = CR.SUdps(icr);
    sortdps = sort(SUdps{:},'descend');
    
    subplot(1,2,1); hold on
    plot(sortdps(1),CR.dprime(icr),'.k')
    
    subplot(1,2,2); hold on
%     plot((sortdps(1)+(mean(sortdps)-sortdps(1)))+(sortdps(2)+(mean(sortdps)-sortdps(2)))+(sortdps(3)+(mean(sortdps)-sortdps(3))),CR.dprime(icr),'.k')
    plot(4*(sum(sortdps)+(mean(sortdps)-sum(sortdps))),CR.dprime(icr),'.k')
end

plot([0 3],[0 3],':k')
axis square

subplot(1,2,1); hold on
plot([0 3],[0 3],':k')
axis square
