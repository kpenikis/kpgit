function compareSoundRanges(SUBJECT,SESSION)
%
%  KP, 2019-12
%

close all

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

% VS
L_VS = sound2dB_Sp(S_VS);

[min(S_VS) max(S_VS)]
[min(L_VS) max(L_VS)]

figure; 
plot(L_AM,'LineWidth',2)
hold on
plot(L_VS,'LineWidth',2)


end