function MC_randxsim
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
%  KP, 2019-01
%

close all

varPar       = 'TrShuff';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'RS'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = {'rand' 'sim'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 500;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
convwin      = exp(-lambda*(1:500));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubpopStart  = 1;

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
q = load(fullfile(rootdir,whichStim,'Full','each','CR_each.mat'));
CReach = q.CR;
clear q

% Check that matching data files were imported
if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
    keyboard
end
if size(Cell_Time_Trial_Stim,1)<size(CReach,1)
    keyboard
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

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];

% Set figsavedir
figsavedir = fullfile(rootdir,whichStim,varPar,'SessionsMin3cells');
if ~exist(fullfile(figsavedir,'backupTables'),'dir')
    mkdir(fullfile(figsavedir,'backupTables'))
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
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>2);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for isess = 1:numel(AllSessions)
    
    if NcellSess(isess)<2
        continue
    end
    whichSess = AllSessions{isess};
    
    % Define cells and stimuli
    [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
        whichSess, nTrialMat, UnitData,...
        theseStim, iRS, iNS, minTrs, convwin, AnWin );
    
    ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);
    
    % Get indices of RS//NS cells (from just "theseCells")
    iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
    iNS = find(UnitInfo(theseCells,:).TroughPeak<0.43 & [UnitData(theseCells).BaseFR]'>2);
    
    
    % Set the subpopulation of cells to use
    UseCells = 1:numel(theseCells);
    SUdps    = CReach.dprime;
    
    switch whichCells
        
        case 'RS'
            UseCells   = iRS;
            SUdps      = CReach(iRS,:).dprime;
        
        case 'BestRS'
            [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
            UseCells   = iRS(iSUdps(1:SubpopStart));
        case 'SkipBestRS'
            [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
            UseCells   = iRS(iSUdps((SubpopStart+1):end));
            
        case 'Mid20RS'
            [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
            UseCells   = iRS(iSUdps(SubpopStart+(0:19)));
        case 'Mid10RS'
            [dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
            UseCells   = iRS(iSUdps(SubpopStart+(0:9)));
            SUdps      = dps(SubpopStart+(0:9));
            
            % Individual sessions
        case 'Mid10RS_Apr02'
            isess = find(strcmp({UnitData(theseCells).Session},'Apr02-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Apr11'
            isess = find(strcmp({UnitData(theseCells).Session},'Apr11-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Apr15'
            % probably were'nt enough trials to run SU d'
            return
        case 'Mid10RS_Oct26'
            % only 8 RS cells
            if SubpopStart>8
                error('not enough cells')
            end
            isess = find(strcmp({UnitData(theseCells).Session},'Oct26-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Mar26'
            isess = find(strcmp({UnitData(theseCells).Session},'Mar26-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Mar28'
            isess = find(strcmp({UnitData(theseCells).Session},'Mar28-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Mar30'
            isess = find(strcmp({UnitData(theseCells).Session},'Mar30-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Jan17'
            isess = find(strcmp({UnitData(theseCells).Session},'Jan17-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Jan21'
            isess = find(strcmp({UnitData(theseCells).Session},'Jan21-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
        case 'Mid10RS_Jan25'
            % only 8 RS cells
            if SubpopStart>8
                error('not enough cells')
            end
            isess = find(strcmp({UnitData(theseCells).Session},'Jan25-AM'));
            idxC  = intersect(isess',iRS);
            [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend');
            UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
            SUdps    = dps(find(dps<1,1,'first')+(0:9))
            
        case 'LoudestRS'
            nSpkRS     = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
            [~,iSUFR]  = sort(nSpkRS,'descend');
            UseCells   = iRS(iSUFR(1:SubpopStart));
        case 'SkipLoudestRS'
            nSpkRS     = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
            [~,iSUFR]  = sort(nSpkRS,'descend');
            UseCells   = iRS(iSUFR((SubpopStart+1):end));
        case 'LowFFRS'
            keyboard
            [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
            UseCells   = iRS(iSUdps((SubpopStart+1):end));
    end
    
    
    if length(UseCells)<2
        continue
    end
    
    
    for TrM = 1:numel(PickTrials)
        
        TrMethod = PickTrials{TrM};
        
        % Train and test classifier
        [S_AssMat_NC,~,~] = runSVMclass_SnE( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
            BootstrapN, nStim, Dur, length(UseCells), TrMethod, TrainSize, TestSize, KernelType ); 
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %##########################################################################
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        %% Plot NEURAL results
        
        ConfMat = mean(S_AssMat_NC,3);
        
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
        savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s_%s',varPar,WinEnds-WinBeg+1,TrainSize,tau,whichStim,whichSess,TrMethod);
        
        print(hf,fullfile(figsavedir,savename),'-dpdf')
        
        
        %% Save results to master table
        
        mastertablesavename = sprintf('CR_v%s',varPar);
        thistablesavename   = savename;
        
        CR1 = table;
        CR1.figname  = {savename};
        CR1.Stim     = {whichStim};
        CR1.Cells    = {whichSess};
        CR1.iC       = length(UseCells);
        CR1.Results  = {S_AssMat_NC};
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
        CR1.SUids    = {UseCells};
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
            save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
            
        else % Save updated table
            
            % Concatenate new data
            CR = [CR; CR1];
            save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
            save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
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
