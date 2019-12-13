function MasterClass
% MasterClass (all parameters defined at top of file)
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality. 
%
%
%  KP, 2019-12
%

% close all

whichIrr     = 'AC';

whichCells   = 'NS'; 
PickTrials   = 'rand';

Un_Index     = [40 115];

BootstrapN   = 1000;
KernelType   = 'linear';

tau          = 5;
lambda       = 1/tau;
winlen       = 500;
convwin      = exp(-lambda*(1:winlen));
% convwin      = convwin./convwin(1);

PSTHsize     = 'N-2';
TrainSize    = 1;
TestSize     = 1;
minTrs       = 15 + TrainSize + TestSize;

rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'Cell_Time_Trial_Stim_simtrs'));
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
    Un_Indices = 1:257;
end

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


%% Prepare input data

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Define cells and stimuli
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ACstim  = 1:8;
DBstim  = [1:6 9:11];

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>3);
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[CTTS,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
                                whichIrr, whichCells, nTrialMat, UnitData,...
                                ACstim, DBstim, iRS, iNS, minTrs, convwin );
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Train and test classifier 
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AssMat = runSVMclass( CTTS, BootstrapN, nStim, Dur, nUns, ...
                  PickTrials,TrainSize, TestSize, KernelType );
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%%
% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

% scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
% shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];


ConfMat = mean(AssMat,3);
muPC = mean(diag(ConfMat))*100;
dprime = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1); 


% Plot
hf = figure;
imagesc(ConfMat)
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s/%i/%i trs  |  %s (N=%i)',muPC,dprime,KernelType,whichIrr,PSTHsize,TrainSize,TestSize,whichCells,nUns))
% title(sprintf('%s SVM  (iUn %s, cat)  PC=%0.1f%%  %i/%i trs, conv %i',KernelType,num2str(Un_Index),mean(diag(mean(AssignmentMat,3)))*100,TrainSize,TestSize,length(convwin)))

% Save figure
savedir = fullfile(fn.figs,'StimClassRcorr');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('Results_Train%i_conv%i_%s_%s_%s_RepPredictors',TrainSize,tau,whichIrr,whichCells,PickTrials);
% savename = sprintf('Results_%s_Train%i_conv%i_iUns%icat%i',KernelType,TrainSize,length(convwin),Un_Index);

print(fullfile(savedir,savename),'-dpdf')


%% Save results to master table 

tablesavename = sprintf('CR_%s_%s',whichIrr,whichCells);

CR1 = table;
CR1.figname  = {savename};
CR1.Stim     = {whichIrr};
CR1.Cells    = {whichCells};
CR1.iC       = 0;
CR1.Results  = {AssMat};
CR1.PC       = muPC;
CR1.dprime   = dprime;
CR1.trials   = {PickTrials};
CR1.rcorrOpt = {'norm'};
CR1.nTemp    = {PSTHsize};
CR1.nTrain   = TrainSize;
CR1.nTest    = TestSize;
CR1.conv     = {'exp'};
CR1.tau      = tau;

% Load saved table
clear q;
try
q=load(fullfile(savedir,tablesavename));
end


if ~exist('q','var')
    
    % Save new table
    CR = CR1;
    save(fullfile(savedir,tablesavename),'CR','-v7.3')
    
else % Save updated table
    
    % Concatenate new data
    CR = q.CR;
    CR = [CR; CR1];
    
    save(fullfile(savedir,tablesavename),'CR','-v7.3')
end

keyboard


end
