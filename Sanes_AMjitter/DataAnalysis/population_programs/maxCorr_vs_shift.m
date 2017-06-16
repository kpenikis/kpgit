function maxCorr_vs_shift(subject)

fn = set_paths_directories;
T = readtable(fullfile(fn.processed,['UnitDataSummary_' subject '.txt']));

[~, plotOptions] = setOptions;
set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultTextInterpreter','none')
set(0,'defaultaxesfontsize',plotOptions.labelSize)


plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );

hf1=figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-500 500],'k','LineWidth',0.5)

hf2 = figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-1 1],'k','LineWidth',0.5)

for ij = unique(T.jitter)'
    if strcmp(ij,'placeholder'), continue, end
    
    fT = T( strcmp(T.jitter,ij) & T.depth~=0 ,:);
    
    figure(hf1)
    plot(fT.maxCorr,fT.shftCorr,'o','LineWidth',2,'MarkerSize',10,'Color',ALLcolors(strcmp(strtok(ij,'_'),strtok(plotOptions.colSelect,'_')),:))
    
    figure(hf2)
    plot(fT.maxCorr,fT.VS,'o','LineWidth',2,'MarkerSize',10,'Color',ALLcolors(strcmp(strtok(ij,'_'),strtok(plotOptions.colSelect,'_')),:))
    
end

figure(hf1)
xlim([-0.15 0.15])
ylim([-200 200])
xlabel('maximum average correlation between single trial spike train and envelope')
ylabel(sprintf('best shift (ms)\nnegative values = increasing delay of spikes relative to sound'))
title('maxCorr vs best shift (ms) (only significantly correlated datapoints)')
set(gcf,'PaperOrientation','landscape');
print(hf1,fullfile(fn.processed,'^Population','Corr_vs_Shift'),'-dpdf','-bestfit');


[r,p]=corrcoef(T(T.depth~=0,:).maxCorr,T(T.depth~=0,:).VS,'rows','complete');

figure(hf2)
text(0.1,0.95,sprintf('r = %2.4f\np = %0.3g',r(1,2),p(1,2)),'FontSize',12)
xlim([-0.15 0.15])
ylim([0 1])
xlabel('average correlation at best shift')
ylabel(sprintf('VS'))
title('maxCorr vs VS (only significantly correlated datapoints)')
set(gcf,'PaperOrientation','landscape');
print(hf2,fullfile(fn.processed,'^Population','Corr_vs_VS'),'-dpdf','-bestfit');



end
