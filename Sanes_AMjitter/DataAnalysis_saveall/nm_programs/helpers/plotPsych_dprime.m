function [x_fit,fitted_yes,fitted_dprime,threshold,slope,linehandle] = ...
    plotPsych_dprime(results,dprimemat,options,plotOptions,zFA,indVar)
%[x,fitted_yes,fitted_dprime,threshold,slope] = ...
%    plotPsych_dprime(results,dprimemat,options,plotOptions,zFA)
%
%This function plots psychometric fits after transforming to dprime space.
%Input variables:
%   results: results structure created by psignifit 4
%   dprimemat: M x 2 matrix arranged as [stimulus (dB re: 100%), dprime]
%   options: structure created by setOptions function for
%            generating fits with psignifit
%   plotOptions: structure created by setOptions function for setting plot
%             features
%   zFA: scored FA rate
%
%
%Written by ML Caras Dec 5 2016
% edited by KP, 2017-02





%Calculate threshold and slope
try
    threshold = getThreshold(results,options.threshPC,false); %scaled threshold
    slope = getSlopePC(results,options.threshPC,false); %scaled slope
    if (strcmp(indVar,'depth')&&threshold>0) || (strcmp(indVar,'jitter')&&threshold<0)
        threshold = nan;
    end
catch
    threshold = nan;
    slope = nan;
end



%Establish x values
xlength = max(results.data(:,1))-min(results.data(:,1));
xLow = min(results.data(:,1))- plotOptions.extrapolLength*xlength;
xHigh = max(results.data(:,1))+ plotOptions.extrapolLength*xlength;
x_fit  = linspace(xLow,xHigh,1000);

%Define percent "yes" responses (y values) for fit
fitted_yes = (1-results.Fit(3)-results.Fit(4))*arrayfun(@(x)...
    results.options.sigmoidHandle(x,results.Fit(1),...
    results.Fit(2)),x_fit)+results.Fit(4);

%Define dprime (y values) for fit
fitted_dprime  = sqrt(2)*erfinv(2*fitted_yes-1)- zFA;


%Plot raw dprime values
plot(dprimemat(:,1),dprimemat(:,2),'.',...
    'MarkerSize',plotOptions.dataSize,...
    'Color',plotOptions.dataColor);
hold on;

%Plot dprime fit
linehandle = plot(x_fit,fitted_dprime,'Color', plotOptions.lineColor,...
    'LineWidth',plotOptions.lineWidth);


%Plot threshold lines
ylimits = get(gca,'ylim');
ymin = ylimits(1);
x = [threshold,threshold];
y = [ymin,options.dprimeThresh];
plot(x,y,'-','Color',plotOptions.lineColor);

xlimits = get(gca,'xlim');
xmin = xlimits(1);
x = [xmin,threshold];
y = [options.dprimeThresh,options.dprimeThresh];
plot(x,y,'-','Color',plotOptions.lineColor);


%Format plot
axis tight
set(gca,'FontSize',plotOptions.fontSize)

ylabel('d''','FontName',plotOptions.fontName,...
    'FontSize', plotOptions.labelSize);
xlabel(plotOptions.xLabel,'FontName',plotOptions.fontName,...
    'FontSize', plotOptions.labelSize);

ylim([floor(10*min([0; dprimemat(:,2)]))/10 ceil(max(dprimemat(:,2)))]);
set(gca,'TickDir','out')
box off

end