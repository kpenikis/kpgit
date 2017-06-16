
subject = 'WWWf_244303';
session = 'CM';

fn = set_paths_directories;

datadir = fullfile(fn.processed,subject);
dataname = sprintf('%s_sess-%s_Phys',subject,session);

% Phys_mat = matfile(fullfile(datadir,dataname));

load(fullfile(datadir,sprintf('%s_sess-%s_Info',subject,session)));
Info.sound_rows

datadir = fullfile(fn.processed,subject);
dataname = sprintf('%s_sess-%s_SoundData',subject,session);
load(fullfile(datadir,dataname))


ds=5;
t_win_sound = 1519500:ds:(1519600+87500);

hf=figure; hold on
scrsz = get(0,'ScreenSize');
set(hf,'Position',[1 scrsz(4)/4 scrsz(3) scrsz(4)/4])
plot(log2(SoundData(1,t_win_sound)./max(SoundData(1,t_win_sound))),'LineWidth',1.5)
plot(SoundData(2,t_win_sound)./max(SoundData(2,t_win_sound)).*1.5 - 6.5,'k')
ylim([-9 0.5])
xlim([1 range(t_win_sound)/ds])

savedir = fn.stim;
savename = 'AMstream_rate_rms_ex';

set(hf,'PaperOrientation','landscape');
print(hf,fullfile(savedir,savename),'-dpdf','-bestfit');


% plot(Phys_mat.Phys(9,5.*t_win_sound)./range(Phys_mat.Phys(9,5.*t_win_sound)).*2 - 2)





