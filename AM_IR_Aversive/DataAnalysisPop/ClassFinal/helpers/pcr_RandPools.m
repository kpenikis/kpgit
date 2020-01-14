

figure;
set(gcf,'Position',widesmall)
hold on
for icr = 1:size(CR,1)
    SUdps = CR.SUdps(icr);
    plotSpread(SUdps{:},'distributionIdx', CR.iC(icr)*ones(size(SUdps{:})),'showMM',3)
end

plot(CR.iC,CR.dprime,'k-')


figure; 
set(gcf,'Position',widesmall)
hold on
for icr = 1:size(CR,1)
    SUdps = CR.SUdps(icr);
    
    subplot(1,2,1); hold on
    plot(max(SUdps{:}),CR.dprime(icr),'ok')
    
    subplot(1,2,2); hold on
    plot(sum(SUdps{:}),CR.dprime(icr),'ok')
end

plot([0 3],[0 3],':k')
axis square

subplot(1,2,1); hold on
plot([0 3],[0 3],':k')
axis square
