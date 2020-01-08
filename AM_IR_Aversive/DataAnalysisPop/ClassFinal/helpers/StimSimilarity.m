
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
winlen       = 500;
convwin      = exp(-lambda*(1:winlen));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AnWin = 501:1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMStim  = 1:8;
SpStim  = 1:8;


%%
% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
smallscreen = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
tallscreen  = [1 scrsz(4) scrsz(3)/4 scrsz(4)];


%% Load data

fn = set_paths_directories('','',1);

% AM DATA

savedir_AM = fullfile(fn.figs,'ClassAM');
rawdata = 'CTTS_AM';
% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir_AM,'RawData',rawdata)); %Cell_Time_Trial_Stim
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% FILTER AM
ETTS_AM = Env_Time_Trial_Stim(:,AnWin,:,AMStim);
clear q Env_Time_Trial_Stim


% SPEECH DATA

savedir = fullfile(fn.figs,'ClassSpeech');
rawdata = 'CTTS_Speech_nonSim';
% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'RawData',rawdata)); %Cell_Time_Trial_Stim
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% FILTER SPEECH
ETTS_Sp = Env_Time_Trial_Stim(:,AnWin,:,SpStim);
clear q Env_Time_Trial_Stim


%% Plot each stimulus

% AM
yspace = 0.8;
hfAM = figure;
set(hfAM,'Position',tallscreen)
hold on
for is = 1:size(ETTS_AM,4)
    plot((is-1)*yspace+mean(mean(ETTS_AM(:,:,:,is),3,'omitnan'),1,'omitnan'),...
        'Color',0.03*[1 1 1],'LineWidth',4)
end
set(gca,'ytick',0:yspace:(yspace*7),'yticklabel',1:8,'xlim',[0 500],'ylim',[0 yspace*8])

print_eps_kp(hfAM,fullfile(savedir_AM,'StimEnvs'))


% Speech
yspace = 0.025;
hfSp = figure;
set(hfSp,'Position',tallscreen)
hold on
for is = 1:size(ETTS_Sp,4)
    plot((is-1)*yspace+mean(mean(ETTS_Sp(:,:,:,is),3,'omitnan'),1,'omitnan'),...
        'Color',0.03*[1 1 1],'LineWidth',4)
end
set(gca,'ytick',0:yspace:(yspace*7),'yticklabel',1:8,'xlim',[0 500],'ylim',[0 yspace*8])

print_eps_kp(hfSp,fullfile(savedir,'StimEnvs'))


%% Correlations between stimuli
keyboard


EnvCorrs_AM = nan(numel(AMStim));
for ist1 = 1:numel(AMStim)
    for ist2 = (ist1+1):numel(AMStim)
        EnvCorrs_AM(ist1,ist2) = corr(mean(mean(ETTS_AM(:,:,:,ist1),3,'omitnan'),1,'omitnan')',...
            mean(mean(ETTS_AM(:,:,:,ist2),3,'omitnan'),1,'omitnan')');
    end
end

EnvCorrs_Sp = nan(numel(SpStim));
for ist1 = 1:numel(SpStim)
    for ist2 = (ist1+1):numel(SpStim)
        EnvCorrs_Sp(ist1,ist2) = corr(mean(mean(ETTS_Sp(:,:,:,ist1),3,'omitnan'),1,'omitnan')',...
            mean(mean(ETTS_Sp(:,:,:,ist2),3,'omitnan'),1,'omitnan')');
    end
end


figure; 
set(gcf,'Position',fullscreen)

subplot(1,3,1)
imagesc(EnvCorrs_AM)
cmocean('algae')
set(gca,'clim',[-0.5 1])
axis square
colorbar
title('AM')

subplot(1,3,2)
imagesc(EnvCorrs_Sp)
cmocean('algae')
set(gca,'clim',[-0.5 1])
axis square
colorbar
title('Speech')

subplot(1,3,3)
hold on
histogram(EnvCorrs_AM,'BinEdges',linspace(-0.5,1,30),'FaceColor','k')
histogram(EnvCorrs_Sp,linspace(-0.5,1,30),'FaceColor','g')
axis square
set(gca,'Color','none')
title('AM (black) vs Speech (green)')

suptitle('Correlations between stimuli')


% set(gcf,'PaperPosition',
orient(gcf,'landscape')
print(gcf,fullfile(savedir,'StimCorrs'),'-dpdf','-bestfit')


%% Now get similarity estimate for each stimulus


figure; hold on

SimScore_AM = nan(numel(AMStim),1);
for ist = 1:8
    if sum(~isnan([EnvCorrs_AM(ist,:) EnvCorrs_AM(:,ist)']))~=7
        keyboard
    end
    plot(ist,[EnvCorrs_AM(ist,:) EnvCorrs_AM(:,ist)'],'ok','MarkerSize',10)
    SimScore_AM(ist) = mean([EnvCorrs_AM(ist,:) EnvCorrs_AM(:,ist)'],'omitnan');
end

SimScore_Sp = nan(numel(SpStim),1);
for ist = 1:8
    if sum(~isnan([EnvCorrs_Sp(ist,:) EnvCorrs_Sp(:,ist)']))~=7
        keyboard
    end
    plot(ist,[EnvCorrs_Sp(ist,:) EnvCorrs_Sp(:,ist)'],'og','MarkerSize',10)
    SimScore_Sp(ist) = mean([EnvCorrs_Sp(ist,:) EnvCorrs_Sp(:,ist)'],'omitnan');
end



% figure; hold on
plot(SimScore_AM,'k','LineWidth',2)
plot(SimScore_Sp,'g','LineWidth',2)
title('Correlation of each stim with others in set')
xlabel('Stimulus')
ylabel('Corr coeff')
set(gca,'Color','none')
print_eps_kp(gcf,fullfile(savedir,'StimCorrs_2'))


% Sorted 
[m_AM,is_AM] = sort(SimScore_AM);
[m_Sp,is_Sp] = sort(SimScore_Sp);


figure;

subplot(2,1,1)
% hold on
% plot([1 8],[0 0],':k')
plot(m_AM,'k','LineWidth',2)
set(gca,'Color','none','xtick',1:8,'xticklabel',is_AM)
ylabel('mean correlation')
xlabel('AM Stim ID')
ylim(0.3*[-1 1])
xlim([0 9])
grid on


subplot(2,1,2)
% hold on
% plot([1 8],[0 0],':k')
plot(m_Sp,'g','LineWidth',2)
set(gca,'Color','none','xtick',1:8,'xticklabel',is_Sp)
ylabel('mean correlation')
xlabel('Speech Stim ID')
ylim(0.3*[-1 1])
xlim([0 9])
grid on

print_eps_kp(gcf,fullfile(savedir,'StimCorrs_3'))






