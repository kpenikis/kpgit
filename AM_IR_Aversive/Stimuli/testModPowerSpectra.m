function testModPowerSpectra(signal,fs)


% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
widescreen = [1 scrsz(4)/4 scrsz(3) scrsz(4)/4];

%%

% signal = ComplexSine';

AMcutoff = 150; 

L = length(signal);

f = fs*(0:(L/2))'/L;

y = fft(signal);

P2 = abs(y/L);

P1 = P2( 1 : L/2+1 );

P1 = P1.*sqrt(f);

P1(2:end-1) = 2*P1(2:end-1);


figure;
set(gcf,'Position',widescreen)

subplot(1,3,1)
plot(f,P1/max(P1),'LineWidth',2)
xlim([0 AMcutoff])


%%

subplot(1,3,2)

plot(f(f<AMcutoff),envelope(P1(f<AMcutoff),20,'rms')/max(envelope(P1(f<AMcutoff),20,'rms')))

% Doesn't look like MS in Ding et al, even with very long speech excerpt.
% What's wrong, or is some pre-processing necessary?? 


%%

subplot(1,3,3)
plotMPS(signal,fs);
xlim([0 AMcutoff])



