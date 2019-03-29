
uiopen('/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AllSpeech.wav',1)
set(0,'DefaultAxesFontSize',10)


%%

data = data/max(data);
absdata = abs(data) - min(abs(data));

%~~~~~~~~~~~~~~~~~~~
% ENVELOPE METHOD 1: find peaks, must be separated by at least X samples

minpeakdist_ms = 12;
minpeakdist_samps = round(minpeakdist_ms/1000 * fs);

ENV_Peaks = envelope( data + (rand(size(data))*0.006-0.003) ,minpeakdist_samps,'peak');
% ENV_Peaks = ENV_Peaks(fs+1:end-fs);
ENV_Peaks(ENV_Peaks<0) = 0;

% idx = find(ENV_Peaks==0);
% ENV_Peaks([1:idx(1) idx(end):end]) = 0;
% ENV_Peaks(1:idx(1)) = 0;

ENV_Peaks   = ENV_Peaks/max(ENV_Peaks);




%~~~~~~~~~~~~~~~~~~~
% ENVELOPE METHOD 2: Hilbert analytic signal, then low pass filtered

hilbtrans = abs(hilbert(data));

% Lowpass
lpFilt = designfilt('lowpassfir','PassbandFrequency',45, ...
         'StopbandFrequency',60,'DesignMethod','kaiserwin','SampleRate',fs);
% fvtool(lpFilt)

ENV_Hilbert = filtfilt(lpFilt,hilbtrans);
% ENV_Hilbert = ENV_Hilbert(fs+1:end-fs);
ENV_Hilbert(ENV_Hilbert<0) = 0;
ENV_Hilbert = ENV_Hilbert/max(ENV_Hilbert);


%~~~~~~~~~~~~~~~~~~~
% ENVELOPE METHOD 3: rms

% data2 = reshape([absdata; -absdata],1,length(absdata)*2);
ENV_RMS = envelope(data,minpeakdist_samps,'rms');
% ENV_RMS = ENV_RMS(fs+1:end-fs);
ENV_RMS(ENV_RMS<0) = 0;
ENV_RMS     = ENV_RMS/max(ENV_RMS);
% ENV_RMS     = ENV_RMS(1:2:end);

% ENV_RMS = envelope(data2,minpeakdist_samps,'rms');
% 
% ENV_RMSfilt = filtfilt(lpFilt,ENV_RMS);
% ENV_RMSfilt = ENV_RMSfilt/max(ENV_RMSfilt);

% data = data(fs+1:end-fs);




%% Plot data and envelopes

L = length(data);
TimeVec = linspace(0,L/fs,L);

% Plot
figure;
plot(TimeVec,data,'k')
% plot(TimeVec,data_filt/max(data_filt),'k')
hold on
% plot(TimeVec,hilbtrans./max(hilbtrans),'LineWidth',2)
plot(TimeVec,ENV_Hilbert,'LineWidth',1.5)
% plot(TimeVec,ENV_Peaks,'LineWidth',1.5)
% plot(TimeVec,ENV_RMS,'LineWidth',1)

xlim([0 max(TimeVec)])
ylim([-1 1.3])


%%
[S,t,f]  = mtspecgramc(data,[minpeakdist_samps*2 minpeakdist_samps/2],params);


%% Plot Modulation Spectra (chronux) 

fmax = 35; 
params.Fs = fs;
params.fpass = [0.1 fmax];

[S_Hilb,f]  = mtspectrumc(ENV_Hilbert,params);
% [S_Peaks]   = mtspectrumc(ENV_Peaks,params);
% [S_RMS,f]     = mtspectrumc(ENV_RMS,params);

S_Hilb  = S_Hilb/sum(S_Hilb(f<=fmax));
% S_Peaks = S_Peaks/sum(S_Peaks(f<=fmax));
% S_RMS   = S_RMS/sum(S_RMS(f<=fmax));

figure; 
plot(f,S_Hilb,'LineWidth',2)
hold on
% plot(f,S_Peaks,'LineWidth',1.5)
% plot(f,S_RMS,'LineWidth',1)
xlim([0 fmax])
xlabel('Frequency (Hz)')



%% Create vocoded sample from any envelope 

envelope_to_Wav(ENV_Hilbert,'WN',fs,'/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AllSpeechVocoded.wav')



%%  Plot Modulation Spectrum (manual) 

f = fs*(0:100:(L/2))'/L;

y_H = fft(ENV_Hilbert);
P2 = abs(y_H/L);
P1_H = P2( 1 : 100 : L/2+1 );
P1_H = P1_H.*sqrt(f);
P1_H(2:end-1) = 2*P1_H(2:end-1);

y_P = fft(ENV_Peaks);
P2 = abs(y_P/L);
P1_P = P2( 1 : 100 : L/2+1 );
P1_P = P1_P.*sqrt(f);
P1_P(2:end-1) = 2*P1_P(2:end-1);


hf_ms = figure; 
subplot(1,2,1); hold on
ip(1)=plot(f,P1_H/max(P1_H),'LineWidth',2);
ip(2)=plot(f,P1_P/max(P1_P),'LineWidth',2);
xlim([0 32])
xlabel('Modulation Frequency (Hz)')
ylabel('Normalized Power')
legend(ip,{'Hilbert' 'Peaks'},'Location','northeast')
hold off

subplot(1,2,2); hold on
plot(f,(P1_P/max(P1_P))-(P1_H/max(P1_H)),'k','LineWidth',2)
plot([0 32],[0 0],'r','LineWidth',2)
xlim([0 32])
xlabel('Modulation Frequency (Hz)')
ylim([-0.5 0.5])
ylabel('MS(env-peaks) - MS(hilbert)')
hold off


%
figure;
plotMPS(ENV_Peaks,fs);
xlim([0 100])


