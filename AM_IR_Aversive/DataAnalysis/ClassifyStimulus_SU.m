function ClassifyStimulus_SU
% ClassifyStimulus_SU
% 
%  Based on RCorr_SU. Changing decision rule to euclidean distance, which
%  is applied in rc_CalcEucDist [was rc_CalcR]. Middle-man function is still
%  run_RCorr. 
% 
%  KP, 2019-10
%


close all
global fn trMin nIterations StimDur exclOnset 

exclOnset   = 0; 
StimDur     = 500;
if exclOnset>0
    StimDur = StimDur-exclOnset;
end
nIterations = 10;
inclSilence = 0;
TempSize    = 15;
app_str   = '_1trTemp';


%% Load Unit data files

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');
if exclOnset>0
    savedir = fullfile(savedir,'exclOnset');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%% Preallocate
if inclSilence
    PCMat = nan(9,9,numel(UnitData));
else
    PCMat = nan(8,8,numel(UnitData));
end


for iUn = 1:numel(UnitData)
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
        
    
    % Load data files
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 || ~exist('TrialData','var')
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
    
    
    % Find all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    if ~isfield(Info,'artifact')
        keyboard
        [all_TDidx,Ntrials,~] = get_clean_trials(TrialData,[],dBSPL,LP,0);
    else
        [all_TDidx,Ntrials,~] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP,0);
    end
    if length(Ntrials)==1
        continue
    end
    
    allStim = unique(TrialData.trID(all_TDidx));
    
    if sum(Ntrials < trMin)==1
        all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<trMin))  = [];
        allStim(Ntrials<trMin)  = [];
        Ntrials(Ntrials<trMin) = [];
    elseif  sum(Ntrials < trMin)>1
        keyboard
    end
    
    
    AllStimStarttimes = cell(1,9);
    
    for istim = allStim'
        
        st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
        
        %%%  Skip ITI stimuli (should be excluded already)
        TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == 0);
        
        % Get timestamps of onsets and offsets
        clear t2 t3 t_win
        t2 = TrialData.onset(TDidx);
        if exclOnset>0
            t2 = t2+exclOnset;
        end
        t3 = TrialData.offset(TDidx);
        t3 = t2 + StimDur;
        if t3(end)>length(Stream_spikes)
            t3 = t3(1:end-1);
            t2 = t2(1:end-1);
        end
        
        % Store starttimes for each trial
        AllStimStarttimes{istim} = t2;
        
    end %istim
    
    
    % Add silence starttimes for silence
    if inclSilence
        AllStimStarttimes{9} = [TrialData.onset(1):StimDur:(TrialData.offset(1)-StimDur)]';
    else
        AllStimStarttimes = AllStimStarttimes(1:8);
    end
    
    
    PCMat(:,:,iUn) = run_RCorr(Stream_spikes,AllStimStarttimes);
    
    
    
end %iUn


% Save data

save(fullfile(savedir,['PCMat' app_str]),'PCMat','-v7.3')


%% Load data

q=load(fullfile(savedir,['PCMat' app_str]));
PCMat  = q.PCMat;

q=load(fullfile(savedir,['PCMat_tsPop_rcorr' app_str '_' num2str(StimDur) 'ms']));
PCMat_Pop = q.PCMat;


%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];


% hist of Mean PC
% hist of Max PC
% hist of PC per stim

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~   Distribution of PC vals for each stimulus   ~~~~~
hf1 = figure;
set(hf1,'Position',fullscreen,'NextPlot','add')
% plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.8 0],'LineWidth',2)
hold on
plot([0 10],[1 1]./size(PCMat_Pop,1),'Color',[0 1 0.8],'LineWidth',2)


PCs=[];
for stid = 1:size(PCMat,1)
    
    for iUn = 1:size(PCMat,3)
        PCs = [PCs; stid PCMat(stid,stid,iUn)];
    end
    
    % Manually make boxplot
    q5  = quantile(PCs(PCs(:,1)==stid,2),0.05);
    q25 = quantile(PCs(PCs(:,1)==stid,2),0.25);
    q75 = quantile(PCs(PCs(:,1)==stid,2),0.75);
    q95 = quantile(PCs(PCs(:,1)==stid,2),0.95);
    
    plot(stid*[1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
    fill(stid+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
    
end

plotSpread(PCs(:,2),'distributionIdx',PCs(:,1),'distributionColors','k','showMM',3)

% Plot population nSpk vector results
plot(diag(PCMat_Pop)','o','Color',[0 1 0.8],'MarkerSize',12,'LineWidth',4)
% plot(diag(PCMat_Vec)','o','Color',[0 1 0.8],'MarkerSize',12,'LineWidth',4)

xlim([0 10])
ylim([0 1])
ylabel('Percent Assigned Correctly')
set(gca,'xtick',1:9,'xticklabel',[Info.stim_ID_key' {'Silence'}],...
    'Color','none','tickdir','out')
title('Stimulus classification, SU data')


%% Save figure

print_eps_kp(hf1,fullfile(savedir,['PC_perStim' app_str]))



end

