
uiopen('/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AllSpeech.wav',1)
set(0,'DefaultAxesFontSize',10)


%%

data = data/max(data);
absdata = abs(data) - min(abs(data));

ENV_type = 'Hilb';

switch ENV_type
    
    case 'Peaks'
        %~~~~~~~~~~~~~~~~~~~
        % ENVELOPE METHOD 1: find peaks, must be separated by at least X samples
        
        minpeakdist_ms = 12;
        minpeakdist_samps = round(minpeakdist_ms/1000 * fs);
        
        ENV = envelope( data + (rand(size(data))*0.006-0.003) ,minpeakdist_samps,'peak');
        ENV(ENV<0) = 0;
        
        ENV   = ENV/max(ENV);
        
        
    case 'Hilb'
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
        ENV = ENV/max(ENV);
        
        
    case 'RMS'
        %~~~~~~~~~~~~~~~~~~~
        % ENVELOPE METHOD 3: rms
        
        ENV = envelope(data,minpeakdist_samps,'rms');
        ENV(ENV<0) = 0;
        ENV = ENV/max(ENV);
        
end


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


%%
[S,t,f]  = mtspecgramc(data,[minpeakdist_samps*2 minpeakdist_samps/2],params);


%% Plot Modulation Spectra (chronux) 

addpath(genpath('/Users/kpenikis/Documents/MATLAB/chronux_2_12'),1)

fmax = 45; 
params.Fs = fs;
params.fpass = [0.1 fmax];

[S,f]  = mtspectrumc(ENV,params);

S  = S/sum(S(f<=fmax));

figure; 
plot(f,S,'LineWidth',2)
xlim([0 fmax])
xlabel('Frequency (Hz)')



%% Create vocoded sample from any envelope 

envelope_to_Wav(ENV,'WN',fs,'/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AllSpeechVocoded.wav')

%% Save amplitude vector 

t1 = 1; %s
t2 = 8; %s
% plot(ENV((t1*fs):(t2*fs)))
% xlim([0 (t2-t1)*fs])

AmpVec = ENV((t1*fs):(t2*fs));
AmpVec = AmpVec-min(AmpVec);
AmpVec = AmpVec/max(AmpVec);

audiowrite('/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AmpVecTest.wav',AmpVec,fs)
save('/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/AmpVecTest.mat','AmpVec','-v7.3')

envelope_to_Wav(AmpVec,'WN',fs,'/Users/kpenikis/Documents/SanesLab/AUDIO samples/SpeechRecordings/VocodedSpeechTest.wav')


%%  Plot Modulation Spectrum (manual) 

f = fs*(0:100:(L/2))'/L;

y_H = fft(ENV);
P2 = abs(y_H/L);
P1_H = P2( 1 : 100 : L/2+1 );
P1_H = P1_H.*sqrt(f);
P1_H(2:end-1) = 2*P1_H(2:end-1);

y_P = fft(ENV);
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
plotMPS(ENV,fs);
xlim([0 100])


