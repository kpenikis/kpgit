
% Create complex sinusoidal function 
fs = 44100;
L = 4e5;
TimeVec = linspace(0,L/fs,L);

ComplexSine = sin(2*pi*2*TimeVec) + sin(2*pi*4*TimeVec) + sin(2*pi*8*TimeVec)...
         + sin(2*pi*16*TimeVec) + sin(2*pi*32*TimeVec) + sin(2*pi*100*TimeVec);%...
%          + sin(2*pi*1000*TimeVec) + sin(2*pi*2000*TimeVec) + sin(2*pi*4000*TimeVec);
ComplexSine = ComplexSine/max(abs(ComplexSine));

rampdur = round(500/1000*fs); %500 ms
n = linspace(-2,2,rampdur);
ramp = (tansig(n) - min(tansig(n)));
ramp = ramp/max(ramp);

ComplexSine(1:rampdur) = ramp.*ComplexSine(1:rampdur);
ComplexSine(end-rampdur+1:end) = fliplr(ramp).*ComplexSine(end-rampdur+1:end);

figure; 
plot(TimeVec,ComplexSine,'k')
ylim([-1.5 1.5])

filename = '/Users/kpenikis/Documents/MATLAB/ModFilter Sound Tools/SpeechTestSentences/ComplexSinusoid.wav';
audiowrite(filename,ComplexSine,fs)
