function makeVStemplates

fs_orig = 24414.0625;
fs_new  = 1000;
[p,q] = rat(fs_new/fs_orig);

fn = set_paths_directories;

j=load(fullfile(fn.stim,'SpeechStim','stim_AsYou_repeat.mat'));
rep_AsYou  = j.buffer;
j=load(fullfile(fn.stim,'SpeechStim','stim_blibBER_repeat.mat'));
rep_ber    = j.buffer;
j=load(fullfile(fn.stim,'SpeechStim','stim_Please_repeat.mat'));
rep_Please = j.buffer;
j=load(fullfile(fn.stim,'SpeechStim','stim_Trees_repeat.mat'));
rep_Trees  = j.buffer;

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/AAB_265054/AAB_265054_sess-Mar30-VS_TrialData.mat')


%~~~~~~~~~~~~~~~~~~
% ~ ~  AS YOU  ~ ~ 
%~~~~~~~~~~~~~~~~~~
nReps_AsYou  = 6;
excess_ms = mod(length(rep_AsYou),nReps_AsYou);
rep_AsYou = rep_AsYou(1:end-excess_ms);
SegDur = length(rep_AsYou)/nReps_AsYou;

Template_AsYou_us = mean(reshape(rep_AsYou',[SegDur nReps_AsYou]),2)';

% Downsample to 1 kHz
Template_AsYou = resample(Template_AsYou_us,p,q);
Template_AsYou(Template_AsYou<0) = 0;

figure; 
plot(rep_AsYou,'k','LineWidth',6)
hold on
plot(repmat(Template_AsYou_us,1,nReps_AsYou),'m','LineWidth',3)
plot(linspace(1,SegDur,length(Template_AsYou)),Template_AsYou,'y','LineWidth',3)
title('AsYou')


%%% Find segments of Sound Stream with high cross correlation
%%% with threshold at 0.95, found ~52 instances of As You 
% figure;
% plot(Template_AsYou./max(Template_AsYou),'k','LineWidth',4)
% hold on
% 
% its   = 1:length(Template_AsYou);
% istep = 50;
% nlags = 100;
% 
% xc=0;
% while all(xc<0.95)
%     its = its+istep;
%     xc = crosscorr(Template_AsYou./max(Template_AsYou),SoundStream(its)./max(SoundStream(its)),nlags);
% end
% [~,im]=max(xc);
% SSseg = SoundStream((its(1)+im-nlags-1)+[1:length(Template_AsYou)]);
% 
% plot(SSseg./max(SSseg),'LineWidth',2)



%~~~~~~~~~~~~~~~~~~
% ~ ~   -ber   ~ ~ 
%~~~~~~~~~~~~~~~~~~
nReps_ber    = 10;
excess_ms = mod(length(rep_ber),nReps_ber);
rep_ber = rep_ber(1:end-excess_ms);
SegDur = length(rep_ber)/nReps_ber;

Template_ber_us = mean(reshape(rep_ber',[SegDur nReps_ber]),2)';

% Downsample to 1 kHz
Template_ber = resample(Template_ber_us,p,q);
Template_ber(Template_ber<0) = 0;

figure; 
plot(rep_ber,'k','LineWidth',6)
hold on
plot(repmat(Template_ber_us,1,nReps_ber),'m','LineWidth',3)
plot(linspace(1,SegDur,length(Template_ber)),Template_ber,'y','LineWidth',3)
title('-ber')


%~~~~~~~~~~~~~~~~~~
% ~ ~  Please  ~ ~ 
%~~~~~~~~~~~~~~~~~~
nReps_Please = 6;
excess_ms = mod(length(rep_Please),nReps_Please);
rep_Please = rep_Please(1:end-excess_ms);
SegDur = length(rep_Please)/nReps_Please;

Template_Please_us = mean(reshape(rep_Please',[SegDur nReps_Please]),2)';

% Downsample to 1 kHz
Template_Please = resample(Template_Please_us,p,q);
Template_Please(Template_Please<0) = 0;

figure; 
plot(rep_Please,'k','LineWidth',6)
hold on
plot(repmat(Template_Please_us,1,nReps_Please),'m','LineWidth',3)
plot(linspace(1,SegDur,length(Template_Please)),Template_Please,'y','LineWidth',3)
title('Please')


%~~~~~~~~~~~~~~~~~~
% ~ ~  Trees   ~ ~ 
%~~~~~~~~~~~~~~~~~~
nReps_Trees = 6;
excess_ms = mod(length(rep_Trees),nReps_Trees);
rep_Trees = rep_Trees(1:end-excess_ms);
SegDur = length(rep_Trees)/nReps_Trees;

Template_Trees_us = mean(reshape(rep_Trees',[SegDur nReps_Trees]),2)';

% Downsample to 1 kHz
Template_Trees = resample(Template_Trees_us,p,q);
Template_Trees(Template_Trees<0) = 0;

figure; 
plot(rep_Trees,'k','LineWidth',6)
hold on
plot(repmat(Template_Trees_us,1,nReps_Trees),'m','LineWidth',3)
plot(linspace(1,SegDur,length(Template_Trees)),Template_Trees,'y','LineWidth',3)
title('Trees')


%% Save all of the templates

save(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'),...
    'Template_AsYou','Template_ber','Template_Please','Template_Trees','-v7.3')


end


