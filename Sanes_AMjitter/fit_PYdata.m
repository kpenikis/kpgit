function [options,results,zFA] = fit_PYdata(data_to_fit,options)
%[options,results] = find_threshPC(data_to_fit,options)
%
%This function determines the percent correct value at which d' ==
%threshold. We use this value as the value at which psignifit 4 should
%calculate threshold and CIs.
%
%Input variables:
%   data_to_fit: M x 3 matrix arranged as [stimulus value, n_yes, n_trials]
%                        (stimulus values are already in dB re: 100% depth)
%
%   options: structure given by setOptions function for fitting
%
%Output variables:
%   options: updated options structure with percent correct value
%   results: structure generated by psignifit 4
%
%Written by ML Caras Dec 5 2016

addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')

n_fa = data_to_fit(1,2);
n_safe = data_to_fit(1,3);
zFA = sqrt(2)*erfinv(2*(n_fa/n_safe)-1);
zThresh = options.dprimeThresh+zFA;
options.threshPC = (erf(zThresh/sqrt(2))+1)/2;

%Fit the data
results = psignifit(data_to_fit,options);

%Remove fields that are much too large to be worth saving
results = rmfield(results,'Posterior');
results = rmfield(results,'weight');



end