function standardPd_PSTHcomparisons( rasters, session, cluname, clulabel,...
                        stimpars, behav, indVar, is )
                    

% Add to comparison plots:
%  - FR in 50 ms preceding
%  - peak rate (instead of average)
%  - add 2nd and 6th periods
%  - total FR pre period


global fn

binsize = 25; 
titlestr = sprintf('%s %s\n%i-%i Hz | %i dBSPL | BehavState: %s | %s',session,cluname,stimpars(1),stimpars(2),stimpars(3),behav,indVar);

savename = sprintf('%s_%s_%i_%i_%i_%i_%i_%s_%s',session,cluname,is,stimpars(1),stimpars(2),stimpars(3),stimpars(4),behav,indVar);
savenameraster = sprintf('%s_%s_%i_Periodic_%i_%i_%i_%i_%s_%s',session,cluname,is,stimpars(1),stimpars(2),stimpars(3),stimpars(4),behav,indVar);

savedir = fullfile(fn.processed,'PSTH_screening');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


rasters = rasters([rasters.AMdepth]==0.75);
if isempty(rasters), return, end

Jitters = [rasters.jitter];
legstr = cell(numel(Jitters),1);
baseFR = calc_baselineFR(rasters);


% Set options for plot
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );
set(0,'DefaultTextInterpreter','none')

nRsp = ceil((numel(Jitters))/2);

scrsz = get(0,'ScreenSize');
if nRsp>1
    figsize = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
    set(0,'DefaultAxesFontSize',16)
else
    figsize = [1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5];
    set(0,'DefaultAxesFontSize',14)
end


% Create figures
hf = figure; 
set(gcf,'Position',figsize,'NextPlot','add');
axis tight
hold on


hc(1) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 100 ms')
ylabel('mean FR standard pd')

hc(2) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 250 ms')
ylabel('mean FR standard pd')

hc(3) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 100 ms')
ylabel('time (ms) of first peak')

hc(4) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 250 ms')
ylabel('time (ms) of first peak')

hc(5) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 50 ms')
ylabel('mean FR standard pd')

hc(6) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR previous 50 ms')
ylabel('time (ms) of first peak')

hc(7) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR response to entire stimulus thus far')
ylabel('mean FR standard pd')

hc(8) = figure;
set(gcf,'Position',[1 scrsz(4)/2.5 scrsz(3)/2 scrsz(4)/2.5],'NextPlot','add');
axis square
hold on
xlabel('FR response to entire stimulus thus far')
ylabel('time (ms) of first peak')



%% First get periodic stimulus data

idxPdc = find(strcmp([rasters.jitter],'0'));
if idxPdc~=1, keyboard, end

legstr{1} = '0';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculate full smoothed psth
plot_periodicPSTH(rasters,idxPdc,binsize,scrsz,plotOptions,titlestr,savedir,savenameraster)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


try
% Prep for other periodic pd spikes
rV = jitter_LUT(rasters(idxPdc).AMrate,char(rasters(idxPdc).jitter));
tV = cumsum( [ ceil(rasters(idxPdc).AMonset+0.75*1000/rV(1)) 1000./rV(2:end) ]);
catch
    keyboard
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear x y
[x,y,prev250FR,prev100FR,prev50FR,totalFR] = standardPd_getspikes(rasters(idxPdc));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

x_p=[]; y_p=[]; p_p=[];
these_pds = [2 6];
for ipd = 1:numel(these_pds)
    
    this_pd = these_pds(ipd);
    
    x_p = [x_p rasters(idxPdc).x( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) ) - tV(this_pd-1) + 1 ];
    y_p = [y_p rasters(idxPdc).y( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) ) ];
    p_p = [p_p repmat(this_pd, size(rasters(idxPdc).x( (rasters(idxPdc).x >= tV(this_pd-1)) & (rasters(idxPdc).x < tV(this_pd)) )) )];
    
    prev250Pdc(ipd) = sum( (rasters(idxPdc).x >= (tV(this_pd-1)-250)) & (rasters(idxPdc).x < tV(this_pd-1)) ) /250*1000/length(rasters(idxPdc).tr_idx);
    prev100Pdc(ipd) = sum( (rasters(idxPdc).x >= (tV(this_pd-1)-100)) & (rasters(idxPdc).x < tV(this_pd-1)) ) /100*1000/length(rasters(idxPdc).tr_idx);
    prev50Pdc(ipd)  = sum( (rasters(idxPdc).x >= (tV(this_pd-1)-50))  & (rasters(idxPdc).x < tV(this_pd-1)) ) /50*1000/length(rasters(idxPdc).tr_idx);
    totalFRpdc(ipd) = sum( (rasters(idxPdc).x >= 0)   &  (rasters(idxPdc).x < tV(this_pd-1)) ) /(tV(this_pd-1))*1000/length(rasters(idxPdc).tr_idx);
    
end

spikes = [];
spikes  = zeros(max(y),250);
spikes2 = zeros(max(y),250);
spikes6 = zeros(max(y),250);

for iy = 1:max(y)
    
    spikes(iy,x(y==iy)) = 1;
    
    spikes2(iy, x_p(y_p==iy & p_p==2) ) = 1;
    spikes6(iy, x_p(y_p==iy & p_p==6) ) = 1;
    
end


%......................................................................

% Calculate smoothed PSTH
FRsmooth_STANDARD = smoothPSTH_v2(spikes);
SDsmooth_STANDARD = smooth_STD(spikes,binsize);

% [FRsmooth_STANDARD,~,~,~,SDbin_STANDARD,SDbinT_STANDARD] = smoothPSTH_v2(spikes);
FRsmooth_2 = smoothPSTH_v2(spikes2);
SDsmooth_2 = smooth_STD(spikes2,binsize);
FRsmooth_6 = smoothPSTH_v2(spikes6);
SDsmooth_6 = smooth_STD(spikes6,binsize);



%% Now plot PSTH of periods 2 and 6 on top of standard

% figure(hf(1))
% subplot(nRsp,2,1); hold on
% title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
% for it = 1:(max(y)+max(y_p(p_p==2))+max(y_p(p_p==6)))
%     plot([0 250],[it it],'Color',[0.5 0.5 0.5],'LineWidth', 0.5*plotOptions.rasterLineWidth)
% end
% plot(x_p(p_p==2),max(y)+y_p(p_p==2),'r+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% plot(x,y,'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% plot(x_p(p_p==6),max(y)+max(y_p(p_p==2))+y_p(p_p==6),'b+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
% ylim([0 max(y)+max(y_p(p_p==2))+max(y_p(p_p==6))+1])
% 
% figure(hf(2))
% subplot(nRsp,2,1); hold on
% title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
% plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
% plot(FRsmooth_2,'r','LineWidth',plotOptions.lineWidth)
% plot(FRsmooth_6,'b','LineWidth',plotOptions.lineWidth)


% Plot PSTH
figure(hf)
subplot(nRsp,2,1); hold on
title(['periodic, ntr=' num2str(numel(rasters(idxPdc).tr_idx))])
fill([1:250 fliplr(1:250)],...
    [FRsmooth_2+SDsmooth_2 fliplr(FRsmooth_2-SDsmooth_2)],...
    'r','FaceAlpha',0.25,'EdgeColor','none')
fill([1:250 fliplr(1:250)],...
    [FRsmooth_6+SDsmooth_6 fliplr(FRsmooth_6-SDsmooth_6)],...
    'b','FaceAlpha',0.25,'EdgeColor','none')
fill([1:250 fliplr(1:250)],...
    [FRsmooth_STANDARD+SDsmooth_STANDARD fliplr(FRsmooth_STANDARD-SDsmooth_STANDARD)],...
    'k','FaceAlpha',0.3,'EdgeColor','none')
plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
plot(FRsmooth_2,'r','LineWidth',plotOptions.lineWidth)
plot(FRsmooth_6,'b','LineWidth',plotOptions.lineWidth)

if 1==(nRsp*2-1)
    set(gca,'xtick',0:50:250); xlabel('Time (ms)'); ylabel('Spikes/sec')
else
    set(gca,'xtick',[],'ytick',[])
end


%% Get time of first peak

[pks,pkts] = findpeaks(FRsmooth_STANDARD,'MinPeakProminence',range(FRsmooth_STANDARD)*0.5);
% plot(pkts,pks,'*g','LineWidth',5)
[pks2,pkts2] = findpeaks(FRsmooth_2,'MinPeakProminence',range(FRsmooth_2)*0.5);
[pks6,pkts6] = findpeaks(FRsmooth_6,'MinPeakProminence',range(FRsmooth_6)*0.5);

if isempty(pkts)
    pkts = 0;
end
if isempty(pkts2)
    pkts2 = 0;
end
if isempty(pkts6)
    pkts6 = 0;
end


%% Plot correlations

% mean FR
figure(hc(1))
hca(1,1) = plot(prev100FR,mean(FRsmooth_STANDARD),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(1,2) = plot(prev100Pdc(1),mean(FRsmooth_2),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(1,3) = plot(prev100Pdc(2),mean(FRsmooth_6),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(2))
hca(2,1) = plot(prev250FR,mean(FRsmooth_STANDARD),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(2,2) = plot(prev250Pdc(1),mean(FRsmooth_2),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(2,3) = plot(prev250Pdc(2),mean(FRsmooth_6),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(5))
hca(5,1) = plot(prev50FR,mean(FRsmooth_STANDARD),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(5,2) = plot(prev50Pdc(1),mean(FRsmooth_2),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(5,3) = plot(prev50Pdc(2),mean(FRsmooth_6),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(7))
hca(7,1) = plot(totalFR,mean(FRsmooth_STANDARD),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(7,2) = plot(totalFRpdc(1),mean(FRsmooth_2),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(7,3) = plot(totalFRpdc(2),mean(FRsmooth_6),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);


% time of first peak FR
figure(hc(3))
hca(3,1) = plot(prev100FR,pkts(1),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(3,2) = plot(prev100Pdc(1),pkts2(1),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(3,3) = plot(prev100Pdc(2),pkts6(1),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(4))
hca(4,1) = plot(prev250FR,pkts(1),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(4,2) = plot(prev250Pdc(1),pkts2(1),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(4,3) = plot(prev250Pdc(2),pkts6(1),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(6))
hca(6,1) = plot(prev50FR,pkts(1),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(6,2) = plot(prev50Pdc(1),pkts2(1),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(6,3) = plot(prev50Pdc(2),pkts6(1),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);

figure(hc(8))
hca(8,1) = plot(totalFR,pkts(1),'ok','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(8,2) = plot(totalFRpdc(1),pkts2(1),'or','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
hca(8,3) = plot(totalFRpdc(2),pkts6(1),'ob','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);


legstr{2} = 'pdc #2';
legstr{3} = 'pdc #6';


%% Now plot jitters on periodic

for ij = 1:numel(Jitters)
    
    % skip periodic
    if ij==idxPdc
        continue
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    clear x_j y_j spikes FRsmooth
    [x_j,y_j,prev250FR,prev100FR,prev50FR,totalFR] = standardPd_getspikes(rasters(ij));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    spikes = zeros(max(y_j),250);
    for iy = 1:max(y_j)
        spikes(iy,x_j(y_j==iy)) = 1;
    end
    
    % Calculate smoothed PSTH
    FRsmooth_Jit = smoothPSTH_v2(spikes);
    SDsmooth_Jit = smooth_STD(spikes,binsize);
    
    % Plot rasters
%     figure(hf(1))
%     subplot(nRsp,2,ij); hold on
%     title([char(rasters(ij).jitter) ', ntr=' num2str(numel(rasters(ij).tr_idx))])
%     for it = 1:(max(y)+max(y_j))
%         plot([0 250],[it it],'k', 'LineWidth', 0.5*plotOptions.rasterLineWidth)
%     end
%     plot(x,y,'k+','MarkerSize',plotOptions.markerSize, 'LineWidth', plotOptions.rasterLineWidth)
%     plot(x_j,y_j+max(y),'+',...
%         'MarkerSize',plotOptions.markerSize,...
%         'LineWidth', plotOptions.rasterLineWidth,...
%         'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
%     ylim([0 max(y_j)+max(y)+1])
%     
%     
%     % Plot PSTH - mean
%     figure(hf(2))
%     subplot(nRsp,2,ij); hold on
%     plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
%     plot(FRsmooth_Jit,'LineWidth',plotOptions.lineWidth,...
%         'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
%     title([char(rasters(ij).jitter) ', ntr=' num2str(numel(rasters(ij).tr_idx))])
    

    % Plot PSTH - with SEM
    figure(hf)
    subplot(nRsp,2,ij); hold on
    
    fill([1:250 fliplr(1:250)],...
        [FRsmooth_STANDARD+SDsmooth_STANDARD fliplr(FRsmooth_STANDARD-SDsmooth_STANDARD)],...
        'k','FaceAlpha',0.3,'EdgeColor','none')
    fill([1:250 fliplr(1:250)],...
        [FRsmooth_Jit+SDsmooth_Jit fliplr(FRsmooth_Jit-SDsmooth_Jit)],...
        ALLcolors(7,:),'FaceAlpha',0.45,'EdgeColor','none')
    
    plot(FRsmooth_STANDARD,'k','LineWidth',plotOptions.lineWidth)
    plot(FRsmooth_Jit,'LineWidth',plotOptions.lineWidth,...
        'Color',ALLcolors(7,:)); %strcmp(strtok(rasters(ij).jitter,'_'),strtok(plotOptions.colSelect,'_')), : ))
    
    rV = jitter_LUT(rasters(ij).AMrate,char(rasters(ij).jitter));
    title([char(rasters(ij).jitter) '  |  prevPd: ' num2str(rV(3)) '  |  ntr=' num2str(numel(rasters(ij).tr_idx))])
    
    if ij==(nRsp*2-1)
        set(gca,'xtick',0:50:250); xlabel('Time (ms)'); ylabel('Spikes/sec')
    else
        set(gca,'xtick',[],'ytick',[])
    end
    
    
    %% Get time of first peak
    [pks,pkts] = findpeaks(FRsmooth_Jit,'MinPeakProminence',range(FRsmooth_Jit)*0.5);
    if isempty(pkts)
        pkts = 0;
    end
    
    %% Plot correlations
    
    % mean FR
    figure(hc(1))
    hca(1,ij+2) = plot(prev100FR,mean(FRsmooth_Jit),'o',...
        'MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    figure(hc(2))
    hca(2,ij+2) = plot(prev250FR,mean(FRsmooth_Jit),'o',...
        'MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    figure(hc(5))
    hca(5,ij+2) = plot(prev50FR,mean(FRsmooth_Jit),'o',...
        'MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    figure(hc(7))
    hca(7,ij+2) = plot(totalFR,mean(FRsmooth_Jit),'o',...
        'MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    
    % time to first peak
    figure(hc(3))
    hca(3,ij+2) = plot(prev100FR,pkts(1),'o','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    
    figure(hc(4))
    hca(4,ij+2) = plot(prev250FR,pkts(1),'o','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    
    figure(hc(6))
    hca(6,ij+2) = plot(prev50FR,pkts(1),'o','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    
    figure(hc(8))
    hca(8,ij+2) = plot(totalFR,pkts(1),'o','MarkerSize',plotOptions.markerSize,'LineWidth',plotOptions.lineWidth);
    
    
    
    legstr{ij+2} = char(rasters(ij).jitter);
    
end  %ij (jitters)


% Set axis limits to be same for all subplots
figure(hf)
hAllAxes = findobj(hf,'type','axes');
ylims = [min(cellfun(@min,get(hAllAxes,'YLim'))) max(cellfun(@max,get(hAllAxes,'YLim')))];
set(hAllAxes,'YLim',ylims)
set(hAllAxes,'XLim',[0 250])

suptitle(titlestr)

set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
% print(hf,fullfile(savedir,savename),'-depsc')


savesuffix = {'meanFR_100ms' 'meanFR_250ms' 'latency_100ms' 'latency_250ms' 'meanFR_50ms' 'latency_50ms' 'meanFR_allFR' 'latency_allFR'};
for ic = 1%:numel(hc)
    
    figure(hc(ic))
    hAllAxes = findobj(hc(ic),'type','axes');
    ylims = [min(get(hAllAxes,'YLim'))-0.5*min(get(hAllAxes,'YLim')) max(get(hAllAxes,'YLim'))+0.25*max(get(hAllAxes,'YLim'))];
    set(hAllAxes,'YLim',ylims)
    xlims = [min(get(hAllAxes,'XLim'))-0.5*min(get(hAllAxes,'XLim')) max(get(hAllAxes,'XLim'))+0.25*max(get(hAllAxes,'XLim'))];
    set(hAllAxes,'XLim',xlims)
    
    legend(hca(ic,:),legstr,'interpreter','none','Location','bestoutside')
    
    
    linfit = polyfit([hca(ic,:).XData], [hca(ic,:).YData], 1);
    fity = polyval(linfit,[xlims(1) [hca(ic,:).XData]  xlims(2)]);
    [r,p] = corrcoef(fity(2:(end-1)),[hca(ic,:).YData]);
    if p(2,1)<0.1
        plot(xlims,[fity(1) fity(end)],'k-','LineWidth',plotOptions.rasterLineWidth)
    else
        plot(xlims,[fity(1) fity(end)],'k:','LineWidth',plotOptions.rasterLineWidth)
    end
    text( 1.15*(xlims(1)+0.2), 0.9*ylims(2), sprintf('r = %3.3f ', r(2,1)))
    
    suptitle(titlestr)
    
    % Save
    if ~exist(fullfile(savedir,'FRhistoryComparison'),'dir')
        mkdir(fullfile(savedir,'FRhistoryComparison'))
    end
    set(gcf,'PaperPositionMode','auto','PaperOrientation','portrait')
    print(hc(ic),fullfile(savedir,'FRhistoryComparison',[savename '_' savesuffix{ic}]),'-depsc')

end


end





