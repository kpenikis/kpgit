function histo_shifts_atMaxCorr(subject)

fn = set_paths_directories;
T = readtable(fullfile(fn.processed,['UnitDataSummary_' subject '.txt']));


[~, plotOptions] = setOptions;
set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultTextInterpreter','none')
set(0,'defaultaxesfontsize',plotOptions.labelSize)

plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );

hf1=figure; hold on
hf2=figure; hold on

for ij = unique(T.jitter)'
    if strcmp(ij,'placeholder'), continue, end
    if (str2double(strtok(ij,'_'))<50) && str2double(strtok(ij,'_'))~=0, continue, end
    
    fT = T( strcmp(T.jitter,ij) & T.depth~=0 ,:);
    
    figure(hf1)
    histogram(fT.shftCorr,-200:5:200,'FaceColor',ALLcolors(strcmp(strtok(ij,'_'),strtok(plotOptions.colSelect,'_')),:),'FaceAlpha',1)
    
    figure(hf2)
    histogram(fT.maxCorr,-0.15:0.005:0.15,'FaceColor',ALLcolors(strcmp(strtok(ij,'_'),strtok(plotOptions.colSelect,'_')),:),'FaceAlpha',1)
    
end

figure(hf1)
xlim([-200 200])
ylabel('N datapoints')
xlabel(sprintf('best shift (ms)\nnegative values = increasing delay of spikes relative to sound'))
title('shift (ms) yielding max correlation, for all significant datapoints')
set(gcf,'PaperOrientation','landscape');
print(hf1,fullfile(fn.processed,'^Population','bestShift_histo'),'-dpdf','-bestfit');

figure(hf2)
xlim([-0.15 0.15])
ylabel('N datapoints')
xlabel('correlation coefficients')
title('average correlation of spike trains to RMS at best shift')
set(gcf,'PaperOrientation','landscape');
print(hf2,fullfile(fn.processed,'^Population','maxCorr_histo'),'-dpdf','-bestfit');


end
