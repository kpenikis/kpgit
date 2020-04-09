function PopSps_CTTS
% 
% PopSps_CTTS
% 

% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT
% sortBy       = 'peakFRquant_Lat'; % 'FRrange'; 
% sortStim     = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHING
convtype     = 'exp';
tau          = 10;
window       = gausswin(tau);
window       = window-min(window);
convwin      = window/sum(window);

% convtype     = 'exp';   
% tau          = 5;
% lambda       = 1/tau;
% % winlen       = 500;
% convwin      = exp(-lambda*(1:500));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minTrs       = 12;


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


%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallhalf    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
widehalf    = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

% Set figsavedir
figsavedir = fullfile(fn.figs,'PopSps');
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
        theseStim  = [4 3 2 1 5 6 7 8]; %1:size(Cell_Time_Trial_Stim,4);
end
% theseStim = theseStim(1:6);

% CellTypes
flagRS = find(UnitInfo.TroughPeak>0.43);
flagNS = find(UnitInfo.TroughPeak<=0.43);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Define cells and stimuli
[CTTS,theseCells,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, AnWin, convtype );

if size(CTTS,1) ~= size(Cell_Time_Trial_Stim,1)
    keyboard
end

% CellTypes
flagRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
flagNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);

ThisData = permute(mean(CTTS(flagRS,:,:,:),3,'omitnan'),[1 2 4 3]);
% ThisData = ThisData - repmat([UnitData(theseCells(flagRS)).BaseFR]'/1000,[1 500 8]);

ETTS = Env_Time_Trial_Stim(:,AnWin,:,theseStim);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% Population sparseness

hf2 = figure;
set(gcf,'Position',widehalf)

for ist = 1:size(ThisData,3)
    
    Env01 = mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1,'omitnan')./max(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1,'omitnan'));
    
    % Cumulative PSTH
    subplot(4,size(ThisData,3),ist)
%     plot([0 500],[0.6135 0.6135],'k')
    plot(Env01./2+0.5,'k','LineWidth',2)
    hold on
    plot(mean(ThisData(:,:,ist),1,'omitnan')./max(mean(ThisData(:,:,ist),1,'omitnan'))./2+0.4,'b','LineWidth',2)
    
    
    % Population sparseness
    subplot(4,size(ThisData,3),size(ThisData,3)+ist)
    plot(Env01./4+0.75,'k','LineWidth',2)
    hold on
    
    % skip nan units
    gu = sum(isnan(ThisData(:,:,ist)),2)==0;
    
    PopSps = nan(1,size(ThisData,2));
    for ims = 1:size(ThisData,2)
        PopSps(1,ims) = calculateSparseness(ThisData(gu,ims,ist));
    end
    plot(PopSps,'m','LineWidth',2)
    
    ylim([0.5 1])
    
    % Stat
%     [r(ist),p(ist)] = corr(mean(ThisData(:,:,ist),1,'omitnan')',PopSps')
%     [r(ist),p(ist)] = corr(diff(mean(ThisData(:,:,ist),1,'omitnan'))',diff(PopSps)');
    [r(ist),p(ist)] = corr(diff(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1,'omitnan'))',diff(PopSps)');
end

r
p
print_eps_kp(hf2,fullfile(figsavedir,['PopSps_' whichStim '_' convtype '_tau' num2str(tau)]))


end
