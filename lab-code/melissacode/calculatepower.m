function  Powermat = calculatepower(Powermat,stimdata,n_trials,stims,...
    which_stim,MF,fs,duration)
%[Powermat,power_vect] = calculatepower(Powermat,stimdata,n_trials,stims,...
%    which_stim,MF,fs,duration)
%
%This function performs an FFT on spike vectors pulled from the stimdata
%input variable, which contains discrete spike times. FFTs are calculated
%for a single data trial at a time.
%The function calls mtspectrumpt.m, which is part of the
%chronux data analysis package. 
%
%Powermat =  [stimulus, ave power, std power, sem power]
%
%ML Caras Dec 2015



%-------------------------------------------------------------------
%Initialize parameters for power analysis

params.tapers = [5 9]; %[TW K] where TW = time-bandwidth product and K =
%the number of tapers to be used (<= 2TW-1). [5 9]
%are the values used by Merri Rosen for her 2010
%analysis.


params.pad = 2;        %Padding for the FFT. -1 corresponds to no padding,
%0 corresponds to the next higher power of 2 and so
%on. This value will not affect the result
%calculation, however, using a value of 1 improves
%the efficiancy of the function and increases the
%number of frequency bins of the result.


params.fpass = [0 10]; %[fmin fmax]
%Frequency band to be used in calculation.

params.Fs = fs;        %Sampling rate

params.err = [1 .05];  %Theoretical errorbars (p = 0.05). For Jacknknife
%errorbars use [2 p]. For no errorbars use [0 p].

params.trialave = 0;   %If 1, average over trials or channels.

dt=1/params.Fs;        %Sampling time

t=0:dt:duration;       %Time grid for prolates for data1 (???)

fscorr = 1;            %If 1, use finite size corrections.

perdur=1/MF;           %Duration of each period (sec)
%-------------------------------------------------------------------

%Create a vector of zeros for each trial that
%the stimulus was presented
power_vect = zeros(n_trials,1);

%Identify the trials which elicited spiking activity
active_trials = unique(stimdata(:,1));

%For each active trial...
for i = 1:numel(active_trials)
    
    %Pull out the spike times
    spikes = stimdata(stimdata(:,1) == active_trials(i),3);
    
    %Calculate the power across frequencies
    [spectra,f,R,Serr] = mtspectrumpt(spikes,params,fscorr,t);
    
    %Find the index value closest to MF
    target = min(abs(f-MF));
    MFidx = find((abs(f-MF) == target));
    
    %Calculate the power at the MF
    MFpower = spectra(MFidx);
    
    %Add the power for that trial into the vector
    power_vect(i) = MFpower;
    
    
end

avePower = mean(power_vect);
stdPower = std(power_vect);
semPower = stdPower/sqrt(n_trials - 1);

Powermat = [Powermat;stims(which_stim),avePower,stdPower,semPower];









end