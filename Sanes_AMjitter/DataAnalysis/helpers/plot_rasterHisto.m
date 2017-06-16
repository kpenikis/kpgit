function plot_rasterHisto(spikes,UnitData)


% Set options for plot
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );
set(0,'DefaultTextInterpreter','none')

scrsz = get(0,'ScreenSize');
figsize = [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];

xlimits = [1 250];

nTrials = size(spikes,1);



% Plot this stimulus
hf = figure; 
set(gcf,'Position',figsize,'NextPlot','add');
hold on

% raster
subplot(2,1,1); hold on
plot([0 0],[-30 30],'k:', 'LineWidth', 2)

% plot([raster(ks).AMonset raster(ks).AMonset],[-30 30],'k:', 'LineWidth', 2)
% plot([raster(ks).stimDur raster(ks).stimDur],[-30 30],'k:', 'LineWidth', 2)

[raster_y,raster_x] = find(spikes);

% long ticks
plot(  raster_x  ,  raster_y  , 'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
for ii=1:nTrials
    plot([xlimits(1) xlimits(2)], [ii ii], 'k', 'LineWidth', plotOptions.rasterLineWidth)
end
axis tight
set(gca, 'XLim', [xlimits(1) xlimits(2)], 'XTick',[], 'YLim', ([0 1+nTrials]))
ylabel('Trials')
set(gca,'FontSize',figFontSize)
hold off


% psth
subplot(2,1,2); hold on
fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*ymaxval + (ymaxval)/2 ,...
    Wcolor,'EdgeColor','none')
plot( t_beg:bin:t_end  , hist_bin , 'k', 'LineWidth', 2)
plot([0 0],[0 ymaxval],'k:', 'LineWidth', 2)
plot([raster(ks).AMonset raster(ks).AMonset],[0 ymaxval],'k:', 'LineWidth', 2)
plot([raster(ks).stimDur raster(ks).stimDur],[0 ymaxval],'k:', 'LineWidth', 2)
set(gca, 'XLim', [t_beg t_end])
xlabel( 'Time (ms)')
ylabel('Spikes/sec')
ylim([0 ymaxval])
set(gca,'FontSize',figFontSize)
hold off

% Set font sizes
set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)

% Add title
t=suptitle(raster(ks).stim_str);
t.FontSize=0.75*figFontSize;







end


