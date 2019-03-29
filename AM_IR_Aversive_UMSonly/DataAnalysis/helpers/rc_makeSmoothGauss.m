function [ts,GW,num_trial_samples] = rc_makeSmoothGauss(sw,Tspan)
% from MLC

%Gaussian smoothing function
tb = 3200; %end points of filter
ts = 1; %sampling time step in ms
t = -tb:ts:tb; %Create time steps (from endpoint to endpoint, in ts steps)


%Gaussian filter:
%
%f(t) = a*exp(-((t-b)^2)/(2(c^2))) + d
%
%a = height of peak (here set to 1)
%b = position of the center of the peak (here set to 0?)
%c = standard deviation of the peak (also called sigma, or Gaussian RMS
%    width, here set to sw in msec)
%d = the value that the function asymptotically approaches far from the
%    peak (here set to 0)
%
%Removing a,b and d from the notation, and substiuting in sw we have:
%f(t) = exp(-(t^2)/(2(sw^2)))
%
%where t is time
%



%Create the gaussian window
numerator = exp(-(t).^2/(2*sw.^2));

%Sum all values to obtain the normalization factor 
denominator = sum(exp(-(t).^2/(2*sw.^2)));

%Create a normalized Gaussian Window whos individual elements sum to 1
GW=numerator/denominator;

%Number of samples in the gaussian
%num_gauss_samples = numel(GW); 

%Number of samples in the trial
num_trial_samples =floor(Tspan/ts);


end