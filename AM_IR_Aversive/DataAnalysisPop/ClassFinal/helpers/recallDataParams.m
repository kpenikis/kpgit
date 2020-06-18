function [CTTS,theseCells,nUns,Dur,nStim,TrainSize,TestSize,UnitData] = recallDataParams(whichStim,whichCells,minTrs,convwin)
% [CTTS,theseCells,nUns,Dur,nStim,TrainSize,TestSize,UnitData] = recallDataParams(whichStim,whichCells,minTrs,convwin)
%
%  KP, 2019-01
%

% close all

varPar       = 'Full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
if nargin<2
    whichCells   = 'each';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tau          = 1;
% lambda       = 1/tau;
% % winlen       = 500;
% convwin      = exp(-lambda*(1:500));
if nargin<4
    convwin      = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
        rawdata = 'CTTS_Speech_nonSim';
        
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
%         theseStim  = 1:size(Cell_Time_Trial_Stim,4);
        theseStim  = [4 3 2 1 5 6 7 8]; 
end

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define cells and stimuli
[CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    whichCells, nTrialMat, UnitData,...
    theseStim, iRS, iNS, minTrs, convwin, AnWin );


end
