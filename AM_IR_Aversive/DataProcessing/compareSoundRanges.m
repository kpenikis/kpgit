function compareSoundRanges(SUBJECT,SESSION)
%
%  KP, 2019-12
%

% close all

%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s-AM_Info'     ,SUBJECT,SESSION); 
load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s-AM_TrialData',SUBJECT,SESSION); 
q=load(fullfile(fn.processed,SUBJECT,filename));
S_AM  = q.SoundStream;
TD_AM = q.TrialData;

filename = sprintf( '%s_sess-%s-VS_TrialData',SUBJECT,SESSION); 
q=load(fullfile(fn.processed,SUBJECT,filename));
S_VS  = q.SoundStream;
TD_VS = q.TrialData;
clear q



%% Convert RMS signals to dB

% AM 
L_AM = sound2dB_AM(S_AM,mode(TD_AM.SPL));

[min(S_AM) max(S_AM)]
[min(L_AM) max(L_AM)]

setFloor = -500;
L_AM(L_AM<-300) = setFloor;
% L_AM(TD_AM.onset(1):end)

% VS
L_VS = sound2dB_Sp(S_VS);

[min(S_VS) max(S_VS)]
[min(L_VS) max(L_VS)]

figure; 
plot(L_AM,'LineWidth',2)
hold on
plot(L_VS,'LineWidth',2)


%% POWER SPECTRUM

% FFT method #1
% [ssf,fws]     = plotMPS(CatEnv,fs);
% [ssf_R,fws_R] = plotMPS(CatPSTH,fs);

fs=1000; 

% FFT method #2
len_AM  = length(L_AM);
f_AM    = fs*(0:(len_AM/2))'/len_AM;
y_AM  = fft(L_AM');
P2_AM = abs(y_AM/len_AM);
P1_AM = P2_AM( 1 : floor(len_AM/2+1) );
P1_AM = P1_AM.*sqrt(f_AM);
P1_AM(2:end-1) = 2*P1_AM(2:end-1);

len_VS    = length(L_VS);
f_VS    = fs*(0:(len_VS/2))'/len_VS;
y_VS  = fft(L_VS');
P2_VS = abs(y_VS/len_VS);
P1_VS = P2_VS( 1 : floor(len_VS/2+1) );
P1_VS = P1_VS.*sqrt(f_VS);
P1_VS(2:end-1) = 2*P1_VS(2:end-1);




% Plot POWER SPECTRUM
figure; 
xlim([0.25 2^9])
set(gca,'xscale','log')
% plot(ssf,fws./max(fws(ssf>0.25&ssf<64)),'LineWidth',2)
hold on
plot(f_AM,P1_AM/max(P1_AM),'LineWidth',2)
plot(f_VS,P1_VS/max(P1_VS),'LineWidth',2)
set(gca,'xtick',[0.5 1 2 4 8 16 32 2^9])
xlabel('Frequency')


end