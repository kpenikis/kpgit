function RCorr_VocSpeech
% RCorr
% 
% 
% KP, 2019-09
%


close all
global fn trMin nIterations StimDur TempSize

exOnShift   = 0; 
StimDur     = 1000;
nIterations = 2000;
TempSize    = 10;
app_str   = '_10trTemp';


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
spkshift = 0;

% Load RMS of stimuli
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/SharedData/shihab_AMdata_2/VSresponses_Clu404.mat');

%% Playing

% Durations = sum(~isnan(AvgEnv),2);
% nSegments = floor(Durations./1000);
% nStimuli  = sum(nSegments);
% 
% [6 10 10 18 6 6]' ./ (Durations./1000);
% 
% % Figure settings
% set(0,'DefaultTextInterpreter','none')
% set(0,'DefaultAxesFontSize',16)
% 
% scrsz = get(0,'ScreenSize');   %[left bottom width height]
% fullscreen   = [1 scrsz(4) scrsz(3) scrsz(4)];
% 
% % Set colors
% colors = [ ...
%       0 150 200 ;...
%       0 200 150;...
%       0   0   0;...
%       0   0   0;...
%     181   0  52;...
%     255  87  51]./255;
% 
% % Print envelopes of speech stimuli
% figure; 
% set(gcf,'Position',fullscreen)
% hold on
% for is = 1:size(AvgEnv,1)
%     plot(is+AvgEnv(is,:),'-','Color',colors(is,:),'LineWidth',3)
% end
% 
% print_eps_kp(gcf,fullfile(fn.figs,'SpeechStimRMS'))


%%

nStim = size(AvgEnv,1);

PCMat = nan(nStim,nStim,numel(UnitData));


for iUn = 1:numel(UnitData)
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
        
    
    % Load data files
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes and shift based on calculated integration time 
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift));  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
        
    end
    
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    %%
    
    [~,~,Stream_spikes,~] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,20,'silence');
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP,0);
    
    allStim = unique(TrialData.trID(all_TDidx));
    
    if sum(Ntrials < trMin*2)==1
        all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<trMin*2))  = [];
        allStim(Ntrials<trMin*2)  = [];
        Ntrials(Ntrials<trMin*2) = [];
    elseif  sum(Ntrials < trMin*2)>1
        keyboard
    end
    
    
    AllStimStarttimes = cell(1,nStim);
    
    for istim = allStim'
        
        % Collect trial indices and timestamps
        
        TDidx = [];
        TDidx = all_TDidx([TrialData.trID(all_TDidx)]==istim);
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        
        if exOnShift
            t2 = t2+250;
        end
        
        % Store starttimes for each trial
        AllStimStarttimes{istim} = t2;
        
    end %istim
    
    PCMat(:,:,iUn) = run_RCorr(Stream_spikes,AllStimStarttimes);
    
    
end %iUn


%% Save data

savedir = fullfile(fn.figs,'RCorr_Speech');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

save(fullfile(savedir,['PCMat_first1000' app_str]),'PCMat','-v7.3')



%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

% hist of Mean PC
% hist of Max PC
% hist of PC per stim

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~   Distribution of PC vals for each stimulus   ~~~~~
hf1 = figure;
set(hf1,'Position',fullscreen,'NextPlot','add')
plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
hold on

PCs=[];
for istim = 1:size(PCMat,1)
    
    for iUn = 1:iUn%size(PCMat,3)
        PCs = [PCs; istim PCMat(istim,istim,iUn)];
    end
    
    % Manually make boxplot
    q5  = quantile(PCs(PCs(:,1)==istim,2),0.05);
    q25 = quantile(PCs(PCs(:,1)==istim,2),0.25);
    q75 = quantile(PCs(PCs(:,1)==istim,2),0.75);
    q95 = quantile(PCs(PCs(:,1)==istim,2),0.95);
    
    plot(istim*[1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
    fill(istim+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
    
end

plotSpread(PCs(:,2),'distributionIdx',PCs(:,1),'distributionColors','k','showMM',3)

xlim([0 nStim+1])
ylim([0 1])
ylabel('Percent Assigned Correctly')
set(gca,'xtick',1:nStim,'xticklabel',[Info.stim_ID_key'],...
    'Color','none','tickdir','out')
title('RCorr classifier, SU data')


%% Save figure

print_eps_kp(hf1,fullfile(savedir,['PC_perStim_first1000' app_str]))



end

