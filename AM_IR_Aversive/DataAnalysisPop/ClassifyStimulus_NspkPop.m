function ClassifyStimulus_NspkPop
% ClassifyStimulus_NspkPop
%
%  Basic linear classifier for segments of Pdc and Irr stimuli.
%
%  Based on RCorr_SU. Changing decision rule to euclidean distance, which
%  is applied in rc_CalcEucDist [was rc_CalcR]. 
%  
%
%  KP, 2019-11 
%


% close all

whichIrr    = 'AC';
TempSize    = 15;
app_str     = '_15trTemp';
nTrials     = 16;

StimDur     = 500;

SetNIterations = 200;

ValidationN = 200;


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in cumulativeSpikeCount)
q=load(fullfile(savedir,'Cell_Time_Trial_Stim_simtrs'));
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
    Un_Indices = 1:257;
end

keyboard
% Find just Apr02 cells
sess_idx = find(strcmp(UnitInfo.Subject,'AAB_265054') & strcmp(UnitInfo.Session,'Apr02-AM'));
Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(intersect(Un_Indices,sess_idx),:,:,:);


%% Preallocate

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end
ACstim  = 1:8;
DBstim  = [1:6 9:11];
ACcells = find(all(nTrialMat(:,ACstim)>=nTrials,2));
DBcells = find(all(nTrialMat(:,DBstim)>=nTrials,2));

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.


%%

switch whichIrr
    case 'AC'
        Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(ACcells,:,:,1:8);
        nTrialMat = nTrialMat(ACcells,ACstim);
    case 'DB'
        Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(DBcells,:,:,[1:6 9:11]);
        nTrialMat = nTrialMat(DBcells,DBstim);
end
nUnits  = size(Cell_Time_Trial_Stim,1);
nStim   = size(Cell_Time_Trial_Stim,4);


%% Bootstrap trials

PCMatN = nan(nStim,nStim,ValidationN);
h = waitbar(0,sprintf('Bootstrapping %i times, grabbing %i trials',ValidationN,nTrials));

for iN = 1:ValidationN
    
    waitbar(iN/ValidationN,h)
    
    % Get nT random trials for each stimulus & cell
    CTTS = nan(nUnits,StimDur,nTrials,nStim);
    for ist = 1:nStim
        for iUn = 1:nUnits
            CTTS(iUn,:,:,ist) = Cell_Time_Trial_Stim(iUn,:,randperm(nTrialMat(iUn,ist),nTrials),ist);
        end
    end
    
    CTTS = Cell_Time_Trial_Stim(:,:,1:nTrials,:);
    
    if TempSize==1 || TempSize==(nTrials-1)
        nIterations = nTrials;
    else
        nIterations = SetNIterations;  %new random draws of trials for templates
    end
    
    
    percentCorrectMat = nan(nStim,nStim,nIterations);
    
    % figure; hold on
    
    for iIt = 1:nIterations
        
        % Choose template and test trials
        if TempSize==1
            temp_trs = iIt;
        elseif TempSize==(nTrials-1)
            temp_trs = 1:nTrials; temp_trs(iIt)=[];
        else
            temp_trs = randperm(nTrials,TempSize);
        end
        test_trs = 1:nTrials; test_trs(temp_trs)=[];
        
        
        % Set templates
        T = cell(2,nStim);
        for ist = 1:nStim
            T(1,ist) = {mean(mean(CTTS(:,:,temp_trs,ist),3),2)'};
        end
        T(2,:) = deal({iIt});
        
        
        % Go through each stimulus
        for isTrue = 1:nStim
            
            % Preallocations
            finalMat = [];
            count = 0;
            
            % Go through each trial
            
            for it = test_trs
                
                S = mean(CTTS(:,:,it,isTrue),2)';
                
                [blockAssignment,maxR] = rc_calcEucDist(T,S); %min Euc dist
                %             [blockAssignment,maxR] = rc_calcR(T,S);      %max rcorr
                count = count+1;
                finalMat(count,:) = [it blockAssignment maxR];
                
            end %it
            
            for isAss =  1:nStim
                percentCorrectMat(isTrue,isAss,iIt) = sum(finalMat(:,2)==isAss)/size(finalMat,1);
            end
            
        end %isTrue
        
    end %iIt
    
    PCMatN(:,:,iN) = mean(percentCorrectMat,3,'omitnan');
    
end %ValidationN

close(h)

PCMat = mean(PCMatN,3,'omitnan');

% Extract statistics from confusion matrix
meanPC = mean(diag(PCMat));
dprime = norminv(mean(diag(PCMat)),0,1) - norminv(mean(PCMat(~diag(diag(PCMat)))),0,1); 


%% Save data

savedir = fullfile(fn.figs,'StimClass');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('PCMat_NspPop_E%s_Apr02',app_str);
save(fullfile(savedir,savename),'PCMatN','-v7.3')


%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

% scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
% shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];

hf2 = figure;
imagesc(PCMat)
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
title(sprintf('%s %s  |  PC=%0.2f, d''=%0.2f',whichIrr,app_str(2:end),meanPC,dprime)) 

% print_eps_kp(hf2,fullfile(savedir,['PCMat_nSpk_Pop' app_str '_' num2str(StimDur) 'ms']))
print(fullfile(savedir,savename),'-dpdf')


end

