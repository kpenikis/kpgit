function MC_bySess
% MasterClass (all parameters defined at top of file)
%  Can process AM or Speech data.
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality.
%
%
%  KP, 2020-01
%

close all

whichClass   = 'Full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'allRS'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = {'rand' 'sim'};% 'rand'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 500;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncells_sp    = 2;

rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        rootdir = fullfile(fn.figs,'ClassAM');
        rawdata = 'CTTS_AM';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassSpeech');
        rawdata = 'CTTS_Speech_sim';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
end

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(rootdir,'RawData',rawdata)); %Cell_Time_Trial_Stim
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% Load SU classification results
q = load(fullfile(rootdir,whichStim,whichClass,'each','CR_each.mat'));
CReach = q.CR;
clear q

% Check that matching data files were imported
if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
    keyboard
end
if size(Cell_Time_Trial_Stim,1)<size(CReach,1)
    keyboard
end


% For rate only classifier: shuffle spiketimes within trial
ShuffOpt = 0;
if strcmp(whichClass,'OnlyRate')
    ShuffOpt = 1;
    convwin  = 1;
end


%% Find sessions with simultaneously recorded units to analyze

[AllSessions,~,CellSessID]= unique({UnitData.Session});
NcellSess = histcounts(CellSessID);


%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

% scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
% tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
% widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];


% Set figsavedir
figsavedir = fullfile(rootdir,whichStim,whichClass,'Sess',whichCells);
if ~exist(figsavedir,'dir')
    mkdir(figsavedir)
end


%% Prepare to parse data

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

switch whichStim
    case 'AC'
        theseStim  = 1:8;
    case 'DB'
        theseStim  = [1:6 9:10];
    case 'Speech'
        theseStim  = 1:size(Cell_Time_Trial_Stim,4);
end


% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define cells and stimuli

[CTTS,idxAllCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,...
    theseStim, iRS, iNS, minTrs, convwin, AnWin );

ETTS = Env_Time_Trial_Stim(idxAllCells,AnWin,:,theseStim);

% Get indices of RS//NS cells (from just "theseCells")
iRS = find(UnitInfo(idxAllCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(idxAllCells,:).TroughPeak<=0.43);

for iss = 1:numel(AllSessions)
    
    if NcellSess(iss)<ncells_sp
        continue
    end
    whichSess = AllSessions{iss};
    fprintf('session %s...\n',whichSess)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set the subpopulation of cells to use
    
    isess = find(strcmp({UnitData(idxAllCells).Session},whichSess));
    
    switch whichCells
        case 'all'
            UseCells  = isess';
        case 'allRS'
            UseCells  = intersect(isess',iRS);
        case 'allNS'
            UseCells  = intersect(isess',iNS);
    end
    if numel(UseCells)<ncells_sp
        continue
    end
    
    SUdps = CReach(UseCells,:).dprime
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for TrM = 1:numel(PickTrials)
        
        if iscell(PickTrials)
            TrMethod = PickTrials{TrM};
        else
            TrMethod = PickTrials;
        end
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % Train and test classifier
        
        %Sum
%         [S_AssMat,~,~] = runSVMclass_SUM( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
%             BootstrapN, nStim, Dur, length(UseCells), TrMethod, TrainSize, TestSize, KernelType, ShuffOpt );
        
        %Cat
        [S_AssMat,~,~] = runSVMclass_notNorm( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
            BootstrapN, nStim, Dur, length(UseCells), TrMethod, TrainSize, TestSize, KernelType, ShuffOpt );
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        %% Plot NEURAL results
        
        ConfMat = mean(S_AssMat,3);
        
        muPC    = mean(diag(ConfMat))*100;
        dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
        
        ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);
        
        % Plot
        hf = figure;
        imagesc(ConfMat)
        axis square
        caxis([-1 1])
        cmocean('curl','pivot',0)
        colorbar
        ylabel('True stim')
        xlabel('Assigned')
        set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
        
        title(sprintf('%0.1f%%, d''=%0.2f\nSVM (%s)  |  %s (N=%i)  |  %s',...
            muPC,dprime,whichStim,whichSess,length(UseCells),TrMethod))
        
        
        % Save figure
        savename = sprintf('Res_%s_%s_%s',whichStim,whichSess,TrMethod);
        
        print(hf,fullfile(figsavedir,savename),'-dpdf')
        
        
        %% Save results to master table
        
        mastertablesavename = sprintf('CR_%s_%s',whichClass,whichCells);
        
        CR1 = table;
        CR1.figname  = {savename};
        CR1.Stim     = {whichStim};
        CR1.Cells    = {whichSess};
        CR1.nC       = numel(UseCells);
        CR1.Results  = {S_AssMat};
        CR1.PC       = muPC;
        CR1.dprime   = dprime;
        CR1.WinBeg   = WinBeg;
        CR1.WinEnd   = WinEnds;
        CR1.trials   = {TrMethod};
        CR1.rcorrOpt = {'norm'};
        CR1.nTemp    = {PSTHsize};
        CR1.nTrain   = TrainSize;
        CR1.nTest    = TestSize;
        CR1.conv     = {'exp'};
        CR1.tau      = tau;
        CR1.UnIDs    = {idxAllCells(UseCells)};
        CR1.CRids    = {UseCells};
        CR1.SUdps    = {SUdps};
%         CR1.TrRes    = {TrialResults};
        
        
        % Load saved table
        clear q;
        if ~exist('CR','var')
        try
            q=load(fullfile(figsavedir,mastertablesavename));
            CR = q.CR;
        end
        end
        
        if ~exist('CR','var')
            
            % Save new table
            CR = CR1;
            save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
            
        else % Save updated table
            
            % Concatenate new data
            CR = [CR; CR1];
            save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
        end
        
    end % TrM
end %Session


keyboard

% Sessions: shuffled vs simultaneous
pcr_ShuffRes

% Plot ranked SU vs subpops
pcr_SUvPools

% Plot min max PC and d' as a function of N cells
pcr_minmaxPC

% Plot PC and d' for each stimulus as a function of N cells
pcr_byStim

% Plot SU dprime vs N spikes
pcr_SUdpFR


end
