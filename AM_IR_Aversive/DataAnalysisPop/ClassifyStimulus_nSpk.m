function ClassifyStimulus_nSpk
% ClassifyStimulus_SU
%
%  Based on RCorr_SU. Changing decision rule to euclidean distance, which
%  is applied in rc_CalcEucDist [was rc_CalcR]. Middle-man function is still
%  run_RCorr.
%
%  KP, 2019-10
%


close all
global fn trMin nIterations StimDur exclOnset TempSize

exclOnset   = 0;
StimDur     = 500;
if exclOnset>0
    StimDur = StimDur-exclOnset;
end
nIterations = 1000;
inclSilence = 0;
TempSize    = 1;
app_str   = '_1trTemp';


%% Load data

fn = set_paths_directories('','',1);

% Load spikes data (created in cumulativeSpikeCount)
load(fullfile(fn.figs,'CmSpCount',['Cell_Time_Trial_Stim_' num2str(20) 'trs']))


nTrials = size(Cell_Time_Trial_Stim,3);
nStim   = size(Cell_Time_Trial_Stim,4);


% Preallocate
PCMat = nan(nStim,nStim,size(Cell_Time_Trial_Stim,1));

h = waitbar(0,'Classifying for each cell...');

for iUn = 1:size(Cell_Time_Trial_Stim,1)
    
    waitbar(iUn/size(Cell_Time_Trial_Stim,1),h)
    
    % continue editing from here
%     ntUsed = zeros(1,nStim);
%     percentAssignment = nan(nStim,nStim,nTrials);
    percentCorrectMat = nan(nStim,nStim,nTrials);
        
    for iIt = 1:nTrials
                
        % Set templates
        T = cell(2,nStim);
        T(1,:) = num2cell(permute(sum(Cell_Time_Trial_Stim(iUn,1:StimDur,iIt,:),2),[1,4,2,3]));
        T(2,:) = deal({iIt});
        
        % Go through each stimulus
        for isTrue = 1:nStim
            
            % Preallocations
            finalMat = [];
            count = 0;
            
            % Go through each trial
            test_trials = 1:nTrials;
            test_trials(iIt) = [];
            for it = test_trials
                
                S = sum(Cell_Time_Trial_Stim(iUn,1:StimDur,it,isTrue),2);
                
                [blockAssignment,maxR] = rc_calcEucDist(T,S); %now minE
                count = count+1;
                finalMat(count,:) = [it blockAssignment maxR];
                
                
            end %it
            
            %percentCorrectMat = nan(nstim,nstim,numel(sws),nIterations);
            for isAss =  1:nStim
                percentCorrectMat(isTrue,isAss,iIt) = sum(finalMat(:,2)==isAss)/size(finalMat,1);
            end
            
        end %isTrue
        
    end %iIt
    
    PCMat(:,:,iUn) = mean(percentCorrectMat,3,'omitnan');
        
end %iUn



%% Save data

savedir = fullfile(fn.figs,'StimClass');
if exclOnset>0
    savedir = fullfile(savedir,'exclOnset');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = ['PCMat_nSpk' app_str '_' num2str(StimDur) 'ms'];
save(fullfile(savedir,savename),'PCMat','-v7.3')



q=load(fullfile(savedir,['PCMat_nSpk' app_str '_' num2str(StimDur) 'ms']));
PCMat = q.PCMat;

% Load population vector classification results
q=load(fullfile(savedir,['PCMat_nSpk_Pop' app_str '_' num2str(StimDur) 'ms']));
PCMat_Pop = q.PCMat;
q=load(fullfile(savedir,['PCMat_nSpk_Vec' app_str '_' num2str(StimDur) 'ms']));
PCMat_Vec = q.PCMat;


%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~   Distribution of PC vals for each stimulus   ~~~~~
hf1 = figure;
set(hf1,'Position',fullscreen,'NextPlot','add')
plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
hold on

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
plot(diag(PCMat_Pop)','o','Color',[1 0.8 0],'MarkerSize',12,'LineWidth',4)
plot(diag(PCMat_Vec)','o','Color',[0 1 0.8],'MarkerSize',12,'LineWidth',4)

xlim([0 7])
ylim([0 1])
ylabel('Percent Assigned Correctly')
set(gca,'xtick',1:size(PCMat,1),...
    'Color','none','tickdir','out')
title(['Stimulus classification from total N spikes over ' num2str(StimDur) 'ms'])


% Save figure

print_eps_kp(hf1,fullfile(savedir,savename))


%%

hf2 = figure;

imagesc(PCMat_Pop)

axis square
caxis([0 1])
cmocean('ice')
colorbar
% print_eps_kp(hf2,fullfile(savedir,['PCMat_nSpk_Pop' app_str '_' num2str(StimDur) 'ms']))
print(fullfile(savedir,['PCMat_nSpk_Pop' app_str '_' num2str(StimDur) 'ms']),'-dpdf')


end

