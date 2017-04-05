function hL = plot_neurometric( dprimes, fitdata, plotOptions, hS, setxlimits, colors )
% hL = plot_neurometric( dprimes, fitdata, plotOptions, hS, setxlimits )
%
%   Plots dprime data and fit. No calculations, just plot.
%   called by:  fit_classifier_data  &  format_formula_data
%
%   KP, 2017-02
%

subplot(hS)
    
hL = nan(numel(dprimes),1);
for ip = 1:numel(dprimes)

    % Plot raw dprime values
    dprimemat = dprimes{ip};
    
    if isempty(dprimemat), continue, end

    plot(dprimemat(:,1),dprimemat(:,2),'.',...
        'MarkerSize',plotOptions.dataSize,...
        'Color',colors(ip,:));
    hold on;
    
    
    % Plot dprime fit
    fitted_dprime = fitdata(ip).fit;
    hL(ip) = plot(fitted_dprime(:,1),fitted_dprime(:,2),'Color', colors(ip,:),...
        'LineWidth',plotOptions.lineWidth);
    if fitdata(ip).p_val > 0.05 
        set(hL(ip),'LineStyle',':');
    end
    
    
    % Format plot
    axis tight
    set(gca,'FontSize',plotOptions.fontSize)
    
    ylabel('d''','FontName',plotOptions.fontName,...
        'FontSize', plotOptions.labelSize);
    xlabel(plotOptions.xLabel,'FontName',plotOptions.fontName,...
        'FontSize', plotOptions.labelSize);
    
    ylim([floor(10*min([0; dprimemat(:,2)]))/10 ceil(max(dprimemat(:,2)))]);
    xlim(setxlimits)
    set(gca,'TickDir','out')
    box off
    
    
    
    % Plot threshold lines
    threshold = fitdata(ip).dprime_threshold;
    ylimits = get(gca,'ylim');
    ymin = ylimits(1);
    x = [threshold,threshold];
    y = [ymin,1];
    plot(x,y,'-','Color',colors(ip,:));
    
    xlimits = get(gca,'xlim');
    xmin = xlimits(1);
    x = [xmin,threshold];
    y = [1 1];
    plot(x,y,'-','Color',colors(ip,:));
    
    
    
    
end
    
end

