function [FRmat,varargout] = calculateFR(stimdata,n,n_trials,FRmat,stims,which_stim)
%[FRmat,spike_vect] = calculateFR(stimdata,n,n_trials,FRmat,stims,which_stim)
%
%Calculates the firing rate (sp/s) for a given stimulus value and appends
%it to a matrix of stimuli and FRs. Firing rates are calculated over a
%given duration of n seconds.
%
%FRmat = [stimulus, aveFR, stdFR, semFR]
%
%ML Caras Dec 2015. Modified Dec 2016 for optional output of spike_vect.


%Create a vector of zeros for each trial that
%the stimulus was presented
spike_vect = zeros(n_trials,1);


%Determine the total number of spikes on each
%trial for which there was spiking activity
numspikes = grpstats(stimdata(:,3),stimdata(:,1),{'numel'});

%Add the spiking activity into the spike vector
spike_vect(1:numel(numspikes)) = numspikes;

%Divide each spike count by the analyzed timescale (n seconds)
spike_vect = spike_vect/n; %sp/s

%Calculate the firing rate
aveFR = mean(spike_vect); %sp/s
stdFR = std(spike_vect);
semFR = stdFR/sqrt(n_trials - 1);

FRmat = [FRmat;stims(which_stim),aveFR,stdFR,semFR];

if nargout>1
    varargout{1} = spike_vect;
end


end