function [finalMat,percentCorrectMat] = Rcorr(SpikeMatrix,template,sw,numReps)
%Written by MLC 11/23/2013.
%Calculates Rcorr between individual spike trains and template spike trains.
%Assigns spike trains to the level that gave a template
%with the highest RCORR.
%Runs through different smoothing window sizes to calculate the optimal
%window size, which gives the highest percent correct.
%SpikeMatrix is an MxN structure containing spike time vectors



%%%%-------------------------------------------------------------
%INITIALIZE PARAMETERS
%%%--------------------------------------------------------------

%Identify song stimulus
wavName = SpikeMatrix(1,1).stim;
if wavName{1}(1) == 'B'
    wavDirectory ='/Users/Melissa/Documents/UW/Manuscripts/J Neurosci_2012/Final Data/Songs/48K_BreedingSongs';
elseif wavName{1}(1) == 'N'
    wavDirectory ='/Users/Melissa/Documents/UW/Manuscripts/J Neurosci_2012/Final Data/Songs/48K_NonBreedingSongs';
end
Wavloadname=[wavDirectory,'/',wavName{1}];
[song,fs,bits] = wavread(Wavloadname);

%Define song duration
Tspan = (length(song)/fs)*1000; %msec
clear song;
clear bits;

%Define delay before stimulus onset
Tdelay = 2000; %msec


%Preallocate percent correct matrix
percentCorrectMat = zeros(numReps,2);

%Create gaussian smoothing window
%ts = (time sample)- length of one sample (in msec)
%GW = Gaussian window normalized to sum to 1, centered at 0, with std = sw
%num_trial_samples = the number of samples in the trial
[ts,GW,num_trial_samples] = makeSmoothGauss_rubelab(sw,Tspan);

%Smooth template spike trains using convolve function
T = smoothTemplates_rubelab_alt(template,Tdelay,ts,num_trial_samples,GW);

%Preallocations
nlevels = size(SpikeMatrix,2);
finalMat = [];
count = 0;


%For each stimulus level...
for i = 1:nlevels
    
    %Identify the number of trials
    ntrials = size(SpikeMatrix(1,i).trial,2);
     
    %Define the level
    level = SpikeMatrix(1,i).level;
    
    %For each trial...
    for n=1:ntrials
        
        %If the trial is a non-template trial...
        if n~=template(1,i).trial
            
            %Smooth spike train
            S = smoothTrains_rubelab_alt(n,i,Tdelay,SpikeMatrix,ts,num_trial_samples,GW);
           
            
            
            %Calculate RCORR values and make level assignment
            [levelAssignment,maxR] = calcR_rubelab(template,S,T);
            count = count+1;
            finalMat(count,:) = [level,n,levelAssignment,maxR];
            
            
        end
    end
end

clear GW;

%Calculate percent correct
correct = find(finalMat(:,1) == finalMat(:,3));
percentCorrect = (numel(correct)/length(finalMat))*100;
percentCorrectMat = [sw,percentCorrect];

end






