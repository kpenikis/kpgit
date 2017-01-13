

x=[0:0.1:2];

pcols = copper(numel(x));

ic=0;
figure; hold on
for ix = x
    ic=ic+1;
    plot([2.^(2-ix) 2.^(2+ix)],[0 0],'.','Color',pcols(ic,:),'MarkerSize',20)
    
end

set(gca,'XScale','log')
plot([0 0],[0 4],'k-')
plot([4 4],[-1 1],'k:')
xlabel('AM rates (Hz)')
xlim([0 16])
ylim([-0.5 0.5])
set(gca,'XTick',[0 1 2 4 8 16])
hold off

set(gcf,'PaperOrientation','landscape');






