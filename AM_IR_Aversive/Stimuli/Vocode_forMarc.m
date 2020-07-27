function Vocode_forMarc(flag_LPfilt)
% Vocode_forMarc
% 
% This function applies single-channel vocoding to any sound. The amplitude
% envelope of the original signal is extracted and the fine-structure is 
% replaced by a random white noise carrier. Input and output data formats
% are wav files.
% 
% Customize these 2 variables:
%  -fs_output -- set as desired, but not too low.
%  -pathname_in -- contains the wav files of sounds to be vocoded.
% 
% The output signals are in the range of -1 to 1. The relative amplitude of
% each sound is preserved by normalizing according to the maximum amplitude
% across all sounds. Post-processing is needed to calibrate the output 
% files to be delivered at a desired SPL.
%
% Optional input variable: flag_LPfilt
% By default, flag_LPfilt is set to 0 (skipping the extra filtering step).
% You probably don't want to apply this filter. The most pure definition of
% a sound envelope signal comes from the Hilbert transform alone. 
% 
% 
% Kristina Penikis, 2020-06
% kpenikis@gmail.com
%

close all


% Set the desired sampling frequency for the output signal
% If you leave this commented, it will default to the original fs of the
% input wav files.
% fs_output = 24414;


% Specify the folder with the original sound recording(s), to be vocoded
pathname_in = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/Stimuli';
files = dir(fullfile(pathname_in,'*.wav'));

pathname_out = fullfile(pathname_in,'Vocoded');
if ~exist(pathname_out,'dir')
    mkdir(pathname_out)
end


%%

if nargin<1
    flag_LPfilt = 0;
end

rng('shuffle')

% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',14)
set(groot,'defaultAxesTickDir', 'out');



%%
% First, step through each wav file and extract envelope signal
% in order to normalize all signals into valid range while retaining the
% relative amplitudes across sounds

maxima=[];

for ii = 1:numel(files)
    
    [data,fs_orig] = audioread(fullfile(pathname_in,files(ii).name));
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Envelope Extraction: Hilbert analytic signal
    ENV = abs(hilbert(data));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Optional lowpass filter (you probably don't need this)
    if flag_LPfilt
        
        lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
            'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',fs_orig);
        %         fvtool(lpFilt)
        
        ENV = filtfilt(lpFilt,ENV);
        
    end
    
    ENV(ENV<0) = 0;
    
    
    % Save maximum of this envelope
    maxima = [maxima max(ENV)];
    
end 



%% 

% Make sure fs_output is set
if ~exist('fs_output','var')
    fs_output = fs_orig;
end


% Step through again; this time vocode and save the output signals

for ii = 1:numel(files)
    
    [data,fs_orig] = audioread(fullfile(pathname_in,files(ii).name));
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Envelope Extraction: Hilbert analytic signal
    ENV = abs(hilbert(data));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Optional lowpass filter
    if flag_LPfilt
        
        lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
            'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',fs_orig);
        %         fvtool(lpFilt)
        
        ENV = filtfilt(lpFilt,ENV);
        
    end
    
    ENV(ENV<0) = 0;
    
    % Normalize envelopes so the max value across all of them equals 1
    ENV = ENV/max(maxima);
    
    
    % Plot data and envelopes
    L = length(data);
    TimeVec = linspace(0,L/fs_orig,L);
    
    figure;
    plot(TimeVec,data,'k')
    hold on
    plot(TimeVec,ENV,'LineWidth',2)
    xlim([0 max(TimeVec)])
    ylim([-1 1.3])
    
    
    % Resample 
    TimeVec_new = linspace(0,max(TimeVec),max(TimeVec)*fs_output);
    ENV_new     = interp1( TimeVec, ENV, TimeVec_new, 'linear');
    
    
    % Vocode
    CarrierNoise  = 2*rand(size(ENV_new))-1;
    
    VocodedSignal = CarrierNoise .* ENV_new;
    
    
    % Add to plot
    plot(TimeVec_new,VocodedSignal)
    
    
    % Save output signal
    savename = sprintf('%s_voc',strtok(files(ii).name,'.'));
    % mat
%     save(fullfile(pathname_out,[savename '.mat']),'VocodedSignal','-v7.3')
    % wav
    audiowrite(fullfile(pathname_out,[savename '.wav']),VocodedSignal,fs_output)
    
end 


end


