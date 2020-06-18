clear; close all; 
addpath('Stimuli')


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',14)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');   %[left bottom width height]
smallfig  = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];


%========
Nch = 8;
%========


% [x,fs] = audioread('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/VocodedSpeechStimuli/RoxaneGay_Bravery.wav');
% [x,fs] = audioread('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Stimuli/SpeechRecordings/VocodedSpeechStimuli/BlibberBlubber_cut.wav');
filedir = '/Users/kristinapenikis/Desktop/talk figs/talk sounds/vocding example';
[x,fs_old] = audioread(fullfile(filedir,'stim_Blibber_sentence.wav'));

FS = 16e3;
x = resample(x, FS, fs_old);

[ vocoded_x,carriers,envelopes] = Vocoder(x, FS,       Nch,      45,      'NOISE',      0);
%[vocoded_x,carriers,envelopes] = Vocoder(x, rate, nchan, cutoff , vocoder_type, verbose)
sound(vocoded_x,FS)

newcolors = cmocean('thermal',size(envelopes,2)+1);

hf=figure;
set(gcf,'Position',smallfig)
hold on
set(gca,'ColorOrder',newcolors)
plot((1:length(x))/FS, envelopes,'LineWidth',2);

xlim([0 2.5])
ylim([0 max(max(envelopes))])
set(gca,'tickdir','out','ticklength',[0.025 0.025],'ytick',[0 max(max(envelopes))],'Color','none')
xlabel('Time (seconds)')
title([num2str(Nch) ' bands'])
title('Natural')

% staggered_carriers = zeros(size(carriers,1),1);
% gap_samps = round(0.25*FS);
% for ii=1:size(carriers,2)
%     staggered_carriers = staggered_carriers + [zeros(length(staggered_carriers)-numel(carriers(1+gap_samps*(ii-1):end,ii)),1); carriers(1+gap_samps*(ii-1):end,ii)];
% end
% staggered_carriers = [staggered_carriers; sum(carriers(1:6400,:),2)];
% sound(staggered_carriers,FS)

serial_carriers = [];
dur_samps = round(0.5*FS);
for ii=1:size(carriers,2)
    serial_carriers = [serial_carriers; carriers(1:dur_samps,ii)];
end
% serial_carriers = [serial_carriers; sum(carriers(1:6400,:),2)];
% sound(serial_carriers,FS)


% Save data

savename = sprintf('BlibVoc_%ich',Nch);
% savename = sprintf('BlibVoc_Natural');

print_eps_kp(hf,fullfile(filedir,savename))

audiowrite(fullfile(filedir,[savename '.wav']),vocoded_x,FS)
audiowrite(fullfile(filedir,[savename '_carriers.wav']),serial_carriers,FS)



%%
% 
%     hilbtrans = abs(hilbert(x));
%     
%     % Lowpass
%     lpFilt = designfilt('lowpassiir','FilterOrder',2, ...
%         'PassbandFrequency',45,'PassbandRipple',0.2,'SampleRate',FS);
%     %         fvtool(lpFilt)
%     
%     envelope = filtfilt(lpFilt,hilbtrans);
%     envelope = filtfilt(lpFilt,envelope);  %filter second time
%     envelope(envelope<0) = 0;
%     envelopes = envelope(:);

