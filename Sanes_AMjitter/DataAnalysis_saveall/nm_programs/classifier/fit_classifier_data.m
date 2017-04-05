function [fitdata,hF] = fit_classifier_data( classifier_data, stim, title_str, indVar )

global xlimits

% Add path and set options for psignifit
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
        colors = copper(numel(classifier_data.PYdata.data));
        xLabel = plotOptions.xLabel;
        setxlim_idx = 2;
    case 'jitter'
%         colors = winter(numel(classifier_data.PYdata.data));
        ALLcolors = winter(numel(plotOptions.colSelect));
        for ic=1:size(stim,1)
            colors(ic,:) = ALLcolors(strcmpi(stim{ic,1},plotOptions.colSelect),:);
        end
        xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
        setxlim_idx = 1;
end

fitdata = struct;
% For each condition
for ic = 1:numel(classifier_data.PYdata.data)
    
    
    data = classifier_data.PYdata.data{ic};
    if isempty(data), continue, end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [options,results,zFA] = fit_PYdata(data, options);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Update plot options
    plotOptions.dataColor = colors(ic,:);
    plotOptions.lineColor = colors(ic,:);
    plotOptions.xminval = floor(data(1,1));
    
    
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
        plotPsych_dprime(results, classifier_data.dprime.data{ic},...
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
    ntrials{ic} = classifier_data.PYdata.data{ic}(2:end,3);
    
    % Get char for legend
    legtext{ic} = stim{ic,1};
    if isempty(legtext{ic})
        legtext{ic}='';
    end

end




%% Fit dprime data (same as dp from formua on FRs)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fitdata.dprime = fit_calc_dprime_threshold( classifier_data.dprime.data, ntrials );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotOptions.xLabel=[];
hL3 = plot_neurometric( classifier_data.dprime.data, fitdata.dprime, plotOptions, hS(3), xlimits, colors );
title(hS(3),'Classifier dprime output fit with sigmoid (nlinfit)')

% Finish plot
try
legend(hL1(~cellfun(@isempty,legtext)),legtext(~cellfun(@isempty,legtext)),'Interpreter','none')
legend(hL2(~cellfun(@isempty,legtext)),legtext(~cellfun(@isempty,legtext)),'Interpreter','none')
legend(hL3(~cellfun(@isempty,legtext)),legtext(~cellfun(@isempty,legtext)),'Interpreter','none')
catch
keyboard
end
suptitle(title_str)



end



