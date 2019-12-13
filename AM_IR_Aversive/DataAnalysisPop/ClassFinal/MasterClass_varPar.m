function MasterClass_varPar
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

varPar       = 'AnWin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'all'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBegs      = 501+[-50 -25 500];
WinEnds      = WinBegs+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'rand';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichIrr     = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 10;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
winlen       = 500;
convwin      = exp(-lambda*(1:winlen));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 16;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'Cell_Time_Trial_Stim'));
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
    Un_Indices = 1:257;
end

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


%% Prepare to parse data

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

ACstim  = 1:8;
DBstim  = [1:6 9:11];

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>3);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for ii = 1:numel(WinBegs)
    
    AnWin = WinBegs(ii):WinEnds(ii);
    
    % Define cells and stimuli
    [CTTS,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
                whichIrr, whichCells, nTrialMat, UnitData,...
                ACstim, DBstim, iRS, iNS, minTrs, convwin, AnWin );
    
    
    % Train and test classifier
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
muPC    = mean(diag(ConfMat))*100;
dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1); 


% Plot
hf(ii) = figure;
imagesc(ConfMat)
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')

title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %i/%i trs  |  %s (N=%i)',...
    muPC,dprime,KernelType,whichIrr,TrainSize,TestSize,whichCells,nUns))


% Save figure
savedir = fullfile(fn.figs,'StimClassRcorr');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s',varPar,WinBegs(ii),TrainSize,tau,whichIrr,whichCells);

print(fullfile(savedir,savename),'-dpdf')


%% Save results to master table 

tablesavename = sprintf('CR_v%s_%s_%s',varPar,whichIrr,whichCells);

CR1 = table;
CR1.figname  = {savename};
CR1.Stim     = {whichIrr};
CR1.Cells    = {whichCells};
CR1.iC       = 0;
CR1.Results  = {AssMat};
CR1.PC       = muPC;
CR1.dprime   = dprime;
CR1.WinBeg   = WinBegs(ii);
CR1.WinEnd   = WinEnds(ii);
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


end %vary classification parameter

end
