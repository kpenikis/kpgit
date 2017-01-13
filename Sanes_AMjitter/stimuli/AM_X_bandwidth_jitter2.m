function [signal params] = AM_X_bandwidth_jitter2(AM_Hz, Mdepth, fc, octaves, dur, srate, noise_exp, jitter)
%AM_X_BANDWIDTH    AM_rate times bandwidth interaction
%   [signal] = AM_X_bandwidth_jitter2(AM_Hz, Mdepth, fc, octaves, dur, srate, noise_exp)
%   AM_Hz: AM in Hz
%   Mdepth: modulation depth in % (e.g. 80)
%   fc: centre frequency
%   octaves: bandwidth in octaves Hz (e.g. 2 or .5)
%   dur: duration of sound
%   srate: sampling rate
%   noisexp: noise exponent, where f^-noise_exp (e.g. 0 = white noise, 1 = pink noise, etc.)
%   jitter: AM jitter, yes or no
%   for example:  s = AM_X_bandwidth_jitter2(20,90,800,4,3,44100,1,'jitter');

%%% This version improves the range of jitter to be within +-.25*(1/AM_Hz),
%%% chosen from a normal distribution (see e.g. Boemio et al. (2005)
%%% paper) by clipping the tails of the normal distribution at
%%% +-.25*(1/AM_Hz). However, this is more or less equivalent to decreasing
%%% the std of the normal ditsribution in the first place, i.e. making the
%%% normal distribution narrower.
%%%%%% Tobias Overath, Jan 2010



t = 0:1/srate:dur-1/srate;
if mod(length(t),2) == 0    % make sure signal has even number of samples
else
    t = t(1:end-1);
end

% create noise with specified spectrum exponent
freq = 1:length(t)/2;
mag = freq.^-noise_exp;         % magnitude spectrum with f^noise_exp
mag = [mag fliplr(mag)];        % for positive and negative frequencies
phase = 2*pi*rand(1,length(mag)/2-1);   % random phase spectrum
phase = [0 phase 0 -1*fliplr(phase)];
X = mag.*exp(1i.*phase);         % create spectrum
x = real(ifft(X));              % create noise


% sinusoidal AM modulator
Mdepth = Mdepth/100;
if exist ('jitter')
    jitter_check = 0;
    while jitter_check == 0
        
    AM_jitt = randn(1,100*dur*AM_Hz)*(1/AM_Hz/6)+1/AM_Hz;  % create normal distribution with std = 1/8 and mean = 1/AM_Hz; 100 times more values to choose from to help with selecting
    AM_jitt = AM_jitt - mean(AM_jitt) + 1/AM_Hz;    % make sure distribution has mean of 1/AM_Hz because of sampling error
    AM_jitts = AM_jitt(find(AM_jitt < 1/AM_Hz+.25*(1/AM_Hz)));  % only include AM rates below 1/AM_Hz+.25*(1/AM_Hz)
    AM_jitts = AM_jitts(find(AM_jitts > 1/AM_Hz-.25*(1/AM_Hz)));    % only include AM rates above 1/AM_Hz-.25*(1/AM_Hz)
    AM_select = randperm(length(AM_jitts)); % additional randomization step (not really necessary)
    AM_jitter = AM_jitts(AM_select(1:round(dur*AM_Hz)));
    
    ons = round(AM_jitter*srate);
    onscum = [0 cumsum(ons)];
    % check whether resulting length of stimulus is within specified range
    % (+-100 ms) of specified duration of stimulus
    if onscum(end)+AM_jitter(end)*srate > srate*(dur - 0.01)
        if onscum(end)+AM_jitter(end)*srate < srate*(dur + 0.01)
            AM = zeros(1,length(x));
            for k = 1:length(AM_jitter)
                AM(1,onscum(k)+1:onscum(k+1)) = 1 + Mdepth*sin(2*pi*1/AM_jitter(k)*t(1:ons(k))-pi/2);
            end
            jitter_check = 1;
        end
    end
    end
    
%     AM_jitter = randn(1,dur*AM_Hz)*(1/AM_Hz/4)+1/AM_Hz;
%     AM_jitter = AM_jitter - mean(AM_jitter) + 1/AM_Hz;
    
else
    AM = 1 + Mdepth*sin(2*pi*AM_Hz*t-pi/2);
end
if AM_Hz == 0
    AM = ones(1,length(x));
end
% bandwidth filter original noise signal
if octaves == 0
    bandwidth = fc;
    xf = sin(2*pi*fc*t);
else
    bandwidth = [fc*2^-(octaves/2) fc*2^(octaves/2)];   % determine bandwidth on exponential scale
    xf = freqfilt(bandwidth(1), bandwidth(2), x, srate);
end

if length(AM) > length(xf)
    AM = AM(1:length(xf));
elseif length(AM) < length(xf)
    AM = [AM zeros(1,length(xf)-length(AM))];
end
% combine filtered signal with AM modulator
xf_AM = xf .* AM;

% scale signal to .1 rms
rms = sqrt(mean(xf_AM.^2));
xf_AM = .1*xf_AM/rms;       

% window signal with 20 ms ramps at beginning & end of sound
xf_AM = wind(srate, 200, xf_AM);   
signal = xf_AM;

% if exist('jitter')
%     filename = ['AM_BW_' num2str(AM_Hz) '_' num2str(octaves) '_40'];
% else
%     filename = [num2str(AM_Hz) 'Hz_' num2str(fc) 'fc_' num2str(octaves) 'oct_' num2str(noise_exp)];
% end
% wavwrite(signal, srate, 16, [filename '.wav'])
% 
if exist ('jitter')
    params.AM_jitter = AM_jitter;
    params.fc = fc;
    params.noise_exp = noise_exp;
    params.Mdepth = Mdepth;
    params.bandwidth = bandwidth;
%     save ([filename '.mat'], 'AM_jitter', 'fc', 'noise_exp', 'Mdepth', 'bandwidth')
else
    params.fc = fc;
    params.noise_exp = noise_exp;
    params.Mdepth = Mdepth;
    params.bandwidth = bandwidth;
   
%     save ([filename '.mat'], 'fc', 'noise_exp', 'Mdepth', 'bandwidth')
end


% S = fft(signal);
% fbins = srate*[0:length(signal)/2-1]/length(signal);
% figure
% plot(fbins,abs(S(1:length(S)/2)),'.-')
% sound(signal,srate)