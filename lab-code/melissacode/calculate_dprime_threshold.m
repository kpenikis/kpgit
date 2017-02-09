function [xfit,yfit,dprime_threshold,dprime_threshold_15] = calculate_dprime_threshold(mat,numtrials,stims,s_ind)
%[xfit,yfit,dprime_threshold,dprime_threshold_15] = ...
%   calculate_dprime_threshold(mat,numtrials,stims,s_ind)
%
%This function fits neurometric dprime values with a sigmoidal fit using
%the nlinfit command. Fits are plotted.
%
%Input variables:
%       mat: M x 4 matrix arranged as 
%           [stimulus, ave FR (or power), std FR (or power), sem FR (or power)]
%
%       numtrials: vector containing number of trials for each stimulus value
%
%       stims: vector of stimulus values (dB re: 100%)
%
%       s_ind: index for subplot
%
%Output variables:
%       xfit: vector of x values for sigmoidal fit
%       yfit: vector of y values for signmoidal fit
%       dprime_threshold = AM depth (dB re: 100%) at which d' = 1
%       dprime_threshold_15 = AM depth (dB re: 100%) at which d' = 1.5
%
%
%Written by ML Caras Dec 5 2016



%For purposes of fitting, discard any stimulus values that were presented
%fewer than 5 times
[c,ia,ib] = intersect(stims,mat(:,1));
numtrials = numtrials(ia);
mat = mat(numtrials>=5,:);

%If there were no spikes at all, dprime == NaN. Remove those.
mat = mat(~isnan(mat(:,2)),:);

if ~isempty(mat)
    x = mat(:,1);
    y = mat(:,2);
    
    %Create a sigmoidal function:  f = y0 + a/(1 + exp(-(x - x0)/b))
    %Parameters (p):
    %p(1):  y0 = min
    %p(2):   a = max - min
    %p(3):   b = slope
    %p(4):  x0 = x coordinate at inflection point
    f = @(p,x) p(1) + p(2) ./ (1 + exp(-(x-p(3))/p(4)));
    
    
    %Establish s vector of initial coefficients (beta0)
    beta0 = [0 20 50 5];
    
    %Set the maximum number of iterations to 1000
    options = statset('MaxIter',10000);
    
    %Do the fitting
    [p, ~, ~, ~, mse] = nlinfit(x,y,f,beta0,options); %#ok<NASGU>
    xfit = linspace(x(1),x(end),1000);
    yfit = f(p,xfit);
    yfit_corr = f(p,x);
    
    %Calculate p value
    [~, p_val] = corrcoef(y,yfit_corr);
    if numel(p_val) > 1
        p_val = p_val(1,2);
    end
    
    %Find threshold at dprime == 1
    dprime_threshold = calc_neurometric_threshold(xfit,yfit,p_val,1);
    
    %Find threshold at dprime == 1.5
    dprime_threshold_15 = calc_neurometric_threshold(xfit,yfit,p_val,1.5);
    
    
    %-----------------------------------------------
    % Plot data
    %-----------------------------------------------
    subplot(2,2,s_ind);
    
    %Raw
    hraw = plot(x,y,'k.','markersize',25,'linewidth',2); %#ok<NASGU>
    hold on
    
    %Fit (dashed line if invalid)
    hfit = plot(xfit,yfit,'k-','linewidth',2);
    
    if p_val>0.05
        set(hfit,'linestyle','--');
    end
  
    %Plot threshold lines if there was a threshold
    if ~isnan(dprime_threshold)
        ylimits = get(gca,'ylim');
        ymin = ylimits(1);
        x = [dprime_threshold,dprime_threshold];
        y = [ymin,1];
        plot(x,y,'-','Color',[0.5 0.5 0.5]);
        
        xlimits = get(gca,'xlim');
        xmin = xlimits(1);
        x = [xmin,dprime_threshold];
        y = [1,1];
        plot(x,y,'-','Color',[0.5 0.5 0.5]);
        
    end
    
    
    %Format
    xlabel('AM Depth (dB)')
    ylabel('dprime')
    myformat(gca);
    
    if max(y) < 3.5
        set(gca,'ylim',[0 3.5]);
    end
    
    
    
else
    xfit = [];
    yfit = [];
    dprime_threshold = NaN;
    dprime_threshold_15 = NaN;
end



%-------------------------------------
%OLDER METHOD USED FOR ARO 2014 AND SEATTLE 2015 POSTER
% options = fitoptions('gauss1', 'Lower', [0 -Inf 0],'Upper',[Inf, 0, Inf]);
% f = fit(x,y,'gauss1',options);
%
% eval_x = -40:0.1:-0.1;
% fit_vals = f(eval_x);
%-------------------------------------
end