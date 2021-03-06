function bestFitParams = fit_weibull(x,y)

% Fit data with weibull distribution

% FOR NOW: correct x values to be positive
x = x - min(x);

% Estimate alpha param to begin search
[~,ind50] = min(abs( y - (max(y)-min(y))/2 ));


    % fit using dB values for x, or %? what to set as 0 depth?
initialSearchParams = [x(ind50) 3+4*(rand(1)-0.5)];
% while any(initialSearchParams<0)
%     initialSearchParams = [x(ind50)+4*(rand(1)-0.5) 3+4*(rand(1)-0.5)];
% end

nll = @(params) nll_weibull_fit(params(1),params(2), x, y );
% nll = @(params) nll_lognorm_fit(params(1),params(2), x, y );

[bestFitParams] = fminsearch(nll,initialSearchParams,optimset('MaxFunEvals',500));

aaa=234;

end