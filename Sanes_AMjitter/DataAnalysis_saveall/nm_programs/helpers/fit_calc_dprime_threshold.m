function fitdata = fit_calc_dprime_threshold(dprimemat,ntrials)
%[xfit,yfit,dprime_threshold,dprime_threshold_15] = calculate_dprime_threshold(mat,numtrials)
%   
%   called by:  fit_classifier_data  &  format_formula_data
%
%   This function fits neurometric dprime values with a sigmoidal fit using
%   the nlinfit command. Fits are plotted.
%
%   Input variables:
%       mat: M x 2 matrix arranged as
%           [stimulus dprime]
%
%       ntrials: vector containing number of trials for each stimulus value
%
%   Output variables:
%       xfit: vector of x values for sigmoidal fit
%       yfit: vector of y values for signmoidal fit
%       dprime_threshold = AM depth (dB re: 100%) at which d' = 1
%       dprime_threshold_15 = AM depth (dB re: 100%) at which d' = 1.5
%
%   Adapted by KP from MLC, 2017-02
% 


% Go through each condition
fitdata = struct();
for ic = 1:numel(dprimemat)
    
    if isempty(dprimemat{ic}), continue, end
    
    mat = dprimemat{ic};
    ntr = ntrials{ic}(2:end);
    
    
    % Remove stim with no spikes in the response (dprime=nan) 
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
        dprime_threshold = find_dprime_threshold(xfit,yfit,p_val,1);
        
        %Find threshold at dprime == 1.5
        dprime_threshold_15 = find_dprime_threshold(xfit,yfit,p_val,1.5);
        
        
        % Add relevant datat to output struct
        fitdata(ic).f = f;
        fitdata(ic).fit = [xfit' yfit'];
        fitdata(ic).p_val = p_val;
        fitdata(ic).dprime_threshold = dprime_threshold;
        fitdata(ic).dprime_threshold_15 = dprime_threshold_15;
        
    else
        xfit = [];
        yfit = [];
        dprime_threshold = nan;
        dprime_threshold_15 = nan;
    end
    
    
    
end %ic

%-------------------------------------
%OLDER METHOD USED FOR ARO 2014 AND SEATTLE 2015 POSTER
% options = fitoptions('gauss1', 'Lower', [0 -Inf 0],'Upper',[Inf, 0, Inf]);
% f = fit(x,y,'gauss1',options);
%
% eval_x = -40:0.1:-0.1;
% fit_vals = f(eval_x);
%-------------------------------------
end