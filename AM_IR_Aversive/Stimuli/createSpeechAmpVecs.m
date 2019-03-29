function createSpeechAmpVecs

set(0,'DefaultAxesFontSize',10)

fs_new = 24414.0625;
maxima=[];

pathname = '/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/Stimuli';
files = dir(fullfile(pathname,'*.wav'));

for ii = 1:numel(files)
    
    [data,fs] = audioread(fullfile(pathname,files(ii).name));
    
    %%
    %~~~~~~~~~~~~~~~~~~~
    % ENVELOPE METHOD 2: Hilbert analytic signal, then low pass filtered
    
    hilbtrans = abs(hilbert(data));
    
    % Lowpass
    lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
        'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',fs);
    %         fvtool(lpFilt)
    
    ENV = filtfilt(lpFilt,hilbtrans);
    ENV = filtfilt(lpFilt,ENV);  %filter second time
    ENV(ENV<0) = 0;    
    
    
    %% Save maximum of this envelope
    
    maxima = [maxima max(ENV)];
    
    
end %this wav file



for ii = 1:numel(files)
    
    [data,fs] = audioread(fullfile(pathname,files(ii).name));
    
    %%
    %~~~~~~~~~~~~~~~~~~~
    % ENVELOPE METHOD 2: Hilbert analytic signal, then low pass filtered
    
    hilbtrans = abs(hilbert(data));
    
    % Lowpass
    %         lpFilt = designfilt('lowpassfir','PassbandFrequency',45, ...
    %             'StopbandFrequency',60,'DesignMethod','kaiserwin','SampleRate',fs);
    lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
        'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',fs);
    %         fvtool(lpFilt)
    
    ENV = filtfilt(lpFilt,hilbtrans);
    ENV = filtfilt(lpFilt,ENV);
    ENV(ENV<0) = 0;
    
    % Normalize envelopes so the max value across all of them equals 1
    ENV = ENV/max(maxima);
    
    
    %% Plot data and envelopes
    
    L = length(data);
    TimeVec = linspace(0,L/fs,L);
    
    % Plot
    figure;
    plot(TimeVec,data,'k')
    hold on
    plot(TimeVec,ENV,'LineWidth',1.5)
    
    xlim([0 max(TimeVec)])
    ylim([-1 1.3])
    
    
    %% Plot Modulation Spectra (chronux)
    
    % addpath(genpath('/Users/kpenikis/Documents/MATLAB/chronux_2_12'),1)
    %
    % fmax = 45;
    % params.Fs = fs;
    % params.fpass = [0.1 fmax];
    %
    % [S,f]  = mtspectrumc(ENV,params);
    %
    % S  = S/sum(S(f<=fmax));
    %
    % figure;
    % plot(f,S,'LineWidth',2)
    % xlim([0 fmax])
    % xlabel('Frequency (Hz)')
    
    
    %% Resample and save
    
    TimeVec_new = linspace(0,max(TimeVec),max(TimeVec)*fs_new);
    buffer = interp1( TimeVec, ENV, TimeVec_new, 'linear');
    
    plot(TimeVec_new,buffer,'LineWidth',1.5)
    
    
    % Save amplitude vector
    
    save(fullfile(pathname,[strtok(files(ii).name,'.') '.mat']),'buffer','-v7.3')
    
    
end %this wav file


end


