

crAmt = 0.0001;

% Data settings
% fn = set_paths_directories('','',1);
% savedir = fullfile(fn.figs,'ClassSpeech',whichStim,'Full','SessSubpops');

% Sessions = {...
%     'Mid10RS_Apr02';...
%     'Mid10RS_Apr11';...
%     'Mid10RS_Mar26';...
%     'Mid10RS_Mar28';...
%     'Mid10RS_Mar30';...
%     'Mid10RS_Jan17';...
%     'Mid10RS_Jan21';...
%     };

Sessions = {...
    'Mid10RS_Mar29-VS' ;...
    'Mid10RS_Apr14-VS' ;...
    'Mid10RS_Mar28-VS' ;...
    'Mid10RS_Apr11-VS' ;...
    'Mid10RS_Mar30-VS' ;...
    'Mid10RS_Apr07-VS' ;...
    'Mid10RS_Apr10-VS' ;...
    'Mid10RS_Apr15-VS' ;...
    'Mid10RS_Apr09-VS' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varPar       = 'Full';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
% whichCells   = 'Mid10RS_Mar28'; % BestRS SkipBestRS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'sim';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 1000;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SubpopSize   = [10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        rawdata = 'CTTS_Speech';
        
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
q = load(fullfile(rootdir,whichStim,varPar,'each','CR_each.mat'));
CReach = q.CR;
clear q

% Check that matching data files were imported
if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
    keyboard
end
if size(Cell_Time_Trial_Stim,1)<size(CReach,1)
    keyboard
end


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
% figsavedir = fullfile(rootdir,whichStim,varPar,whichCells);
% if ~exist(fullfile(figsavedir,'backupTables'),'dir')
%     mkdir(fullfile(figsavedir,'backupTables'))
% end


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
for ii = 1:numel(Sessions)
    
    whichCells = Sessions{ii};
    
    for nc = 1:numel(SubpopSize)
        
        nCell = SubpopSize(nc);
        
        % Define cells and stimuli
        [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
            'each', nTrialMat, UnitData,...
            theseStim, iRS, iNS, minTrs, convwin, AnWin );
        
        ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);
        
        % Get indices of RS//NS cells (from just "theseCells")
        iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
        iNS = find(UnitInfo(theseCells,:).TroughPeak<0.43 & [UnitData(theseCells).BaseFR]'>2);
        
        
        % Set the subpopulation of cells to use
        switch whichCells
            
            case 'BestRS'
                [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                UseCells   = iRS(iSUdps(1:nCell));
            case 'SkipBestRS'
                [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                UseCells   = iRS(iSUdps((nCell+1):end));
                
            case 'Mid20RS'
                [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                UseCells   = iRS(iSUdps(nCell+(0:19)));
            
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
                if nCell>8
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
                if nCell>8
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
                UseCells   = iRS(iSUFR(1:nCell));
            case 'SkipLoudestRS'
                nSpkRS     = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
                [~,iSUFR]  = sort(nSpkRS,'descend');
                UseCells   = iRS(iSUFR((nCell+1):end));
            case 'LowFFRS'
                keyboard
                [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                UseCells   = iRS(iSUdps((nCell+1):end));
               
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%   Speech
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'Mid10RS_Mar29-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Mar29-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr14-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr14-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Mar28-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Mar28-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr11-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr11-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Mar30-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Mar30-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr07-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr07-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr10-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr10-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr15-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr15-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
            case 'Mid10RS_Apr09-VS'
                isess = find(strcmp({UnitData(theseCells).Session},'Apr09-VS'));
                idxC  = intersect(isess',iRS);
                [dps,iSUdps] = sort(CReach(idxC,:).dprime,'descend')
%                 UseCells = idxC(iSUdps(find(dps<1,1,'first')+(0:9)));
%                 SUdps    = dps(find(dps<1,1,'first')+(0:9))
        end
        
    end
end
keyboard

