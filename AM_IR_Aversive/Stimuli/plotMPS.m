function plotMPS(w,fs)
%Adapted from matlab file exchange script by: Diego Barrag�n G.
%https://www.mathworks.com/matlabcentral/fileexchange/12240-am-modulation-and-its-spectrum

Ts = 1/fs;

n=floor(log(length(w))/log(2));
N=2^n; 
fw=abs(fft (w(1:N)));
ssf=(-N/2:N/2-1)/(Ts*N);

fws=fftshift(fw);
plot(ssf,fws);
xlabel('Frequency')
ylabel('Magnitude')