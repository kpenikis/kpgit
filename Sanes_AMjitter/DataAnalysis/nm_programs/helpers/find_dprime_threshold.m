function dprime_threshold = calculate_dprime_threshold(xfit,yfit,p_val,dprime_val)
%dprime_threshold = calc_neurometric_threshold(xfit,yfit,p_val,dprime_val)
%
%This function finds the stimulus value at which the neurometric fit
%crosses a dprime of dprimeval.
%
%Input variables:
%   xfit: vector of x values for neurometric fit
%   yfit: vector of y values for neurometric fit (must be same size as x)
%   p_val: p value from pearson's r to determine if fit is valid
%   dprime_val: value at which you want to find threshold
%
%Written by ML Caras Dec 5 2016


%If the fit is not valid, there is no threshold
if isnan(p_val) || p_val > 0.05 
    dprime_threshold = nan;
 
%If the fit is valid, and it crossed dprime_val...
elseif max(yfit) > dprime_val && min(yfit) < dprime_val
    
    indmax = find(yfit == max(yfit));
    indmin = find(yfit == min(yfit));
    
    indmax = indmax(1);
    indmin = indmin(end);
    
    xmax = xfit(indmax);
    xmin = xfit(indmin);
    
    %And if the fit slope was positive
    if xmax > xmin
        
        %Find threshold @ dprime = dprime_val
        target = min(abs(yfit - dprime_val));
        thresh_ind = find(abs(yfit - dprime_val) == target);
        dprime_threshold = xfit(thresh_ind(1));
        
    else
        dprime_threshold = nan;
        
    end
    
%If the fit is valid, but the maximum value is still below dprime_val, 
%there is no threshold
elseif max(yfit) < dprime_val
    
    dprime_threshold = nan;
    
    
    
%If the fit is valid, but the minimum value is above dprime_val, set the
%threshold to the lowest AM depth tested that day
elseif min(yfit) > dprime_val
    
    dprime_threshold = min(xfit);
    
end

end
