function PopResp_CTTS
% 
% PopResp_CTTS
% 

% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT
sortBy       = 'peakFRquant_Lat'; 
sortStim     = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHING
convtype     = 'gauss';
tau          = 20;
window       = gausswin(tau);
window       = window-min(window);
convwin      = window/sum(window);

% convtype     = 'exp';   
% tau          = 5;
% lambda       = 1/tau;
% % winlen       = 500;
% convwin      = exp(-lambda*(1:500));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minTrs       = 5;


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
% q = load(fullfile(rootdir,whichStim,'Full','each','CR_each.mat'));
% CReach = q.CR;
% clear q
% 
% % Check that matching data files were imported
% if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
%     keyboard
% end
% if size(Cell_Time_Trial_Stim,1)<size(CReach,1)
%     keyboard
% end


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
figsavedir = fullfile(fn.figs,'PopResp');
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
[CTTS,~,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'PopResp', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, AnWin, convtype );

if size(CTTS,1) ~= size(Cell_Time_Trial_Stim,1)
    keyboard
end
FR_vec = permute(mean(CTTS,3,'omitnan'),[1 2 4 3]);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% Sort data

clear i_sorted
ThisData = FR_vec;

switch sortBy
    case 'peakFRquant_Lat'
        
        ytickset = 0;
        yticklab = {''};
        
        % RS cells
        [pkRS,ipkRS] = max(ThisData(flagRS,1:250,sortStim),[],2);
        Qbounds      = quantile(pkRS,[0 0.2 0.4 0.6 0.8 1]);
        
        dataRS  = [];
        isortRS = [];
        for iq=1:5
            idx = find(pkRS>=Qbounds(iq) & pkRS<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkRS(idx),'descend');
            isortRS  = [isortRS; idx(idx_sort)];
            dataRS   = [dataRS; pkRS(idx(idx_sort)) ipkRS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab ['RS' num2str(iq)]];
        end
        
        
        % NS cells
        [pkNS,ipkNS] = max(ThisData(flagNS,1:250,sortStim),[],2);
        Qbounds      = quantile(pkNS,[0 1]);
        
        dataNS  = [];
        isortNS = [];
        for iq=1
            idx = find(pkNS>=Qbounds(iq) & pkNS<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkNS(idx),'descend');
            isortNS  = [isortNS; idx(idx_sort)];
            dataNS   = [dataNS; pkNS(idx(idx_sort)) ipkNS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab 'NS'];
        end
        
        sortdata     = [dataRS; dataNS];
        i_sorted     = [flagRS(isortRS); flagNS(isortNS)];
        ytickset     = cumsum(ytickset);
end


%% Plot results

hf1 = figure;
% set(gcf,'Position',tallhalf)
if strcmp(whichStim,'AC')
    set(gcf,'Position',fullscreen)
elseif strcmp(whichStim,'Speech')
    set(gcf,'Position',widehalf)
else
    keyboard
end

for ist = 1:size(ThisData,3)
    
    subplot(1,size(ThisData,3),ist)
    thisRespPlot
    
end

print_eps_kp(hf1,fullfile(figsavedir,['PopResp_' whichStim '_' sortBy]))

set(hf1,'PaperOrientation','landscape')
print(hf1,'-dpdf','-r500','-fillpage', fullfile(figsavedir,['PopResp_' whichStim '_' sortBy]))

% set(hf1,'PaperOrientation','portrait')
% print(hf1,'-dpdf','-r500','-fillpage', fullfile(figsavedir,['PopResp_' whichStim '_' sortBy '_portrait']))

end
