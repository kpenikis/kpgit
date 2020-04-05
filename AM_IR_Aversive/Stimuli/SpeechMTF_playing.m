
uiopen('/Users/kpenikis/Desktop/SpeechRecordings/BlibberBlubber_abridged.wav',1)

%%

% Find peaks
absdata = abs(data) - min(abs(data));
ENV = envelope(absdata,120,'peak');
ENV(ENV<0) = 0;

% Hilbert
hilbtrans = hilbert(data);

% Lowpass
lpFilt = designfilt('lowpassfir','PassbandFrequency',32, ...
         'StopbandFrequency',40,'DesignMethod','kaiserwin','SampleRate',fs);
fvtool(lpFilt)

ENV_LP = filter(lpFilt,absdata);


Amplitude = ENV/max(ENV);
L = length(Amplitude);
TimeVec = linspace(0,L/fs,L);

% Plot
figure;
plot(TimeVec,data/max(data),'k')
hold on

fdelay = grpdelay(lpFilt);
fdelay = round(unique(fdelay));
plot(TimeVec(1:end-fdelay+1),ENV_LP(fdelay:end)/max(ENV_LP),'r','LineWidth',2)

hilbtransLP = abs(hilbert(ENV_LP));
plot(TimeVec(1:end-fdelay+1),hilbtransLP(fdelay:end)/max(hilbtransLP),'b','LineWidth',2)


plot(TimeVec,abs(hilbtrans)/max(abs(hilbtrans)),'LineWidth',3)
plot(TimeVec,Amplitude,'LineWidth',3)
xlim([0 max(TimeVec)])


%%

WN = rand(size(Amplitude))*2-1; %white noise

Wp = [ 300  18000] * 2 / fs;        %cutoff fqs (passband)
Ws = [ 100  20000] * 2 / fs;        %cutoff fqs (stopband)
[N,Wnf] = buttord( Wp, Ws, 3, 20);  %create filter parameters
[B,A] = butter(N,Wnf);              %build filter

WNfilt = filtfilt( B, A, WN );

f=1000;
PureTone = sin(2*pi*f*TimeVec)';

%
carrier = WNfilt;
%

output = carrier.*Amplitude;


%% 

signal = WN;output;ENV;abs(hilbtrans);

L = length(signal);

f = fs*(0:(L/2))'/L;

y = fft(signal);

P2 = abs(y/L);

P1 = P2( 1 : L/2+1 );

P1 = P1.*sqrt(f);

P1(2:end-1) = 2*P1(2:end-1);

% AMcutoff = 32; 
figure;
% plot(f(f<AMcutoff),P1(f<AMcutoff)/max(P1(f<AMcutoff)),'LineWidth',2)
plot(f,P1/max(P1),'LineWidth',2)
xlim([0 10000])


%%

figure;
plot(f(f<AMcutoff),envelope(P1(f<AMcutoff),20,'rms')/max(envelope(P1(f<AMcutoff),20,'rms')))

% Doesn't look like MS in Ding et al, even with very long speech excerpt.
% What's wrong, or is some pre-processing necessary?? 


%%

figure;
plotMPS(ENV,fs);
xlim([0 100])





