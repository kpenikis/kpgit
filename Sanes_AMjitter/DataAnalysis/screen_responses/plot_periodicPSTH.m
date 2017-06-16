function plot_periodicPSTH(rasters,idxPdc,binsize,scrsz,plotOptions,titlestr,savedir,savename)

% Set window for full raster/histo
t_full = [-199 200+rasters(idxPdc).stimDur];
t_full(2) = t_full(2)-mod(t_full(2),binsize);

spikes  = zeros( max(rasters(idxPdc).y), diff(rasters(idxPdc).window_ms) );
for iy = 1:max(rasters(idxPdc).y)
    spikes(iy, rasters(idxPdc).x(rasters(idxPdc).y==iy) - rasters(idxPdc).window_ms(1) +1 ) = 1;
end

FRsmooth = smoothPSTH_v2(spikes);
SDsmooth = smooth_STD(spikes,binsize);

FRsmooth = FRsmooth(:, 1+((t_full(1)-rasters(idxPdc).window_ms(1)) : (-rasters(idxPdc).window_ms(1)+t_full(2))) ); 
SDsmooth = SDsmooth(:, 1+((t_full(1)-rasters(idxPdc).window_ms(1)) : (-rasters(idxPdc).window_ms(1)+t_full(2))) );


% Get rounded ymax val
yvals = [50 100 200 300 400];
yin = find((max(yvals,10+max(FRsmooth+SDsmooth)) - yvals)==0);
ymaxval = yvals(yin(1));


% Get sitmulus waveform
Wave = ap_stimplotting('IIIf_230115',rasters(idxPdc));
Wcolor = [0.85 0.85 0.85];
    
% Create figure
hr = figure;
set(gcf,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],'NextPlot','add');
axis tight
hold on

% raster
subplot(2,1,1); hold on
plot([0 0],[-1 30],'k:', 'LineWidth', 2)
plot([rasters(idxPdc).AMonset rasters(idxPdc).AMonset],[-1 30],'k:', 'LineWidth', 2)
plot([rasters(idxPdc).stimDur rasters(idxPdc).stimDur],[-1 30],'k:', 'LineWidth', 2)
plot(rasters(idxPdc).x,rasters(idxPdc).y,'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
for ii=1:max(rasters(idxPdc).y)
    plot(t_full, [ii ii], 'k', 'LineWidth', plotOptions.rasterLineWidth)
end
ylim([0 max(rasters(idxPdc).y)+1])
xlim(t_full)
set(gca,'xtick',[],'ytick',[0 max(rasters(idxPdc).y)])
ylabel('Trials')

% psth
subplot(2,1,2); hold on
fill([1:length(Wave.y) length(Wave.y):-1:1] /Wave.fs*1000, [Wave.y fliplr(-Wave.y)]/4*ymaxval + (ymaxval)/2 ,...
    Wcolor,'EdgeColor','none')
plot([0 0],[-1 ymaxval],'k:', 'LineWidth', 2)
plot([rasters(idxPdc).AMonset rasters(idxPdc).AMonset],[-1 ymaxval],'k:', 'LineWidth', 2)
plot([rasters(idxPdc).stimDur rasters(idxPdc).stimDur],[-1 ymaxval],'k:', 'LineWidth', 2)

fill([t_full(1):t_full(2) fliplr(t_full(1):t_full(2))],...
    [FRsmooth+SDsmooth fliplr(FRsmooth-SDsmooth)],...
    'k','FaceAlpha',0.3,'EdgeColor','none')
plot(t_full(1):t_full(2),FRsmooth,'k','LineWidth',plotOptions.lineWidth)

xlabel( 'Time (ms)')
ylabel('Spikes/sec')
% set(gca,'FontSize',plotOptions.fontSize)
ylim([0 ymaxval])
xlim(t_full)
hold off


% Add title
suptitle(titlestr)


% Save figure
set(gcf,'PaperPositionMode','auto','PaperOrientation','landscape')
print(hr,fullfile(savedir,savename),'-depsc')



end