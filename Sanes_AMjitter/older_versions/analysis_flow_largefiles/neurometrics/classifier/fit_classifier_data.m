function [fitdata,hF] = fit_classifier_data( classifier_data, title_str, indVar )

global xlimits

% Add path and set options for psignifit
addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')
[options, plotOptions] = setOptions;
options.dprimeThresh = 1;


%% Fit PY data (with process analogous to psychometric process)

% Set up current figure/subplot handles
hF = figure; hold on
scrsz = get(0,'ScreenSize');
set(hF,'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2],...
    'Nextplot','add');
for isp = 1:3
    hS(isp)=subplot(1,3,isp);
    set(hS(isp),'Nextplot','add');
end

switch indVar
    case 'depth'
        colors = copper(numel(classifier_data.PYdata));
        xLabel = plotOptions.xLabel;
        setxlim_idx = 2;
    case 'jitter'
        colors = winter(numel(classifier_data.PYdata));
        xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
        setxlim_idx = 1;
end

% For each condition
for ic = 1:numel(classifier_data.PYdata)
    
    
    data = classifier_data.PYdata{ic};
    if isempty(data), continue, end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [options,results,zFA] = fit_PYdata(data, options);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Update plot options
    plotOptions.xminval = floor(data(1,1));
    plotOptions.dataColor = colors(ic,:);
    plotOptions.lineColor = colors(ic,:);
    
    
    % Plot the percent correct values and fit, and save handles
    plotOptions.xLabel=[];
    subplot(hS(1)); hold on
    [~,hL1(ic)] = plotPsych(results,plotOptions);
    xlimits = get(hS(1),'XLim'); 
    xlimits(setxlim_idx) = 0;
    set(hS(1),'XLim',xlimits)
    title(hS(1),'Classifier PY estimates, psignifit fit')
    
    % Convert to dprime space
    plotOptions.xLabel = xLabel;
    subplot(hS(2)); hold on
    [x,fitted_yes,fitted_dprime,threshold,slope,hL2(ic)] = ...
        plotPsych_dprime(results, classifier_data.dprime{ic},...
                        options,plotOptions,zFA,indVar);
    set(hS(2),'XLim',xlimits)
    title(hS(2),'Classifier PY estimates, psignifit fit transformed to dprime')
    
    
    % Save data to output structure
    fitdata.PYdata(ic).results = results;
    fitdata.PYdata(ic).options = options;
    fitdata.PYdata(ic).x = x;
    fitdata.PYdata(ic).fitted_yes = fitted_yes;
    fitdata.PYdata(ic).fitted_dprime = fitted_dprime;
    fitdata.PYdata(ic).threshold = threshold;
    fitdata.PYdata(ic).slope = slope;
    
    
    % Build ntrials for the fitting below
    ntrials{ic} = classifier_data.PYdata{ic}(2:end,3);
    
    % Get char for legend
    legtext{ic} = classifier_data.stim{ic,1};
end




%% Fit dprime data (same as dp from formua on FRs)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fitdata.dprime = fit_calc_dprime_threshold( classifier_data.dprime, ntrials );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotOptions.xLabel=[];
hL3 = plot_neurometric( classifier_data.dprime, fitdata.dprime, plotOptions, hS(3), xlimits, colors );
title(hS(3),'Classifier dprime output fit with sigmoid (nlinfit)')

% Finish plot
legend(hL1,legtext,'Interpreter','none')
legend(hL2,legtext,'Interpreter','none')
legend(hL3,legtext,'Interpreter','none')

suptitle(title_str)



end



