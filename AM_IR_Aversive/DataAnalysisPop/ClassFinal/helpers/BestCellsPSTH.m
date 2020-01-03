



close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells       =  30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
winlen       = 500;
convwin      = exp(-lambda*(1:winlen));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AnWin = 501:1000;


%%
% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
smallscreen = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];


%% Load data

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        savedir = fullfile(fn.figs,'ClassAM');
        rawdata = 'CTTS_AM';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
    case 'Speech'
        savedir = fullfile(fn.figs,'ClassSpeech');
        rawdata = 'CTTS_Speech_nonSim';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
end

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'RawData',rawdata)); %Cell_Time_Trial_Stim
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% Check that correct UnitData file was imported
if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
    keyboard
end


% Load SU classification results
q = load(fullfile(savedir,whichStim,'Full','each','CR_each.mat'));
CReach = q.CR;
clear q


%% Get "theseCells"

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

switch whichStim
    case 'AC'
        theseStim  = 1:8;
    case 'DB'
        theseStim  = [1:6 9:11];
    case 'Speech'
        theseStim  = 1:size(Cell_Time_Trial_Stim,4);
end

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>2);


% Define cells and stimuli
[CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,...
    theseStim, iRS, iNS, minTrs, convwin, AnWin );
ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);

iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<0.43 & [UnitData(theseCells).BaseFR]'>2);


% Find CTTS indices of best N SUs (RS only!) 
        
[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
UseCells = iRS(iSUdps(1:Ncells));

figure;
set(gcf,'Position',fullscreen)
plot(1:length(dps),dps,'.k','MarkerSize',20)
xlim([0 length(dps)])
ylim([-0.05 2.5])
set(gca,'ytick',0:0.5:2.5,'xtick',[5 10 30 50 100 150])
grid on
title([whichStim ' -- SU dprimes'])
ylabel('d prime')
xlabel('Cell number')
keyboard
print_eps_kp(gcf,fullfile(savedir,'cdf_SUdps'))


%% Plot PSTHs of best N SUs 

figure;
set(gcf,'Position',fullscreen)
for ist = 1:nStim
    subplot(4,2,ist)
    hold on
    plot(mean(CTTS(UseCells,:,:,ist),3,'omitnan')'.*1000,'k','LineWidth',0.5)
    plot(mean(mean(CTTS(UseCells,:,:,ist),3,'omitnan'),1).*1000,'m','LineWidth',4)
    ylim([0 500])
    title(ist)
end
suptitle([whichStim ' - PSTH best ' num2str(Ncells) ' cells'])
keyboard
print_eps_kp(gcf,fullfile(savedir,['PSTH_best' num2str(Ncells)]))


%% Correlations between stimuli
keyboard

if ~exist('nStim','var')
    nStim = 8;
end
EnvCorrs = nan(nStim);
for ist1 = 1:nStim
    for ist2 = (ist1+1):nStim
        EnvCorrs(ist1,ist2) = corr(mean(mean(Env_Time_Trial_Stim(:,AnWin,:,ist1),3,'omitnan'),1,'omitnan')',...
            mean(mean(Env_Time_Trial_Stim(:,AnWin,:,ist2),3,'omitnan'),1,'omitnan')');
    end
end


figure; 
set(gcf,'Position',fullscreen)

subplot(1,3,1)
imagesc(EnvCorrsAM)
cmocean('algae')
set(gca,'clim',[-0.5 1])
axis square
colorbar
title('AM')

subplot(1,3,2)
imagesc(EnvCorrs)
cmocean('algae')
set(gca,'clim',[-0.5 1])
axis square
colorbar
title('Speech')

subplot(1,3,3)
hold on
histogram(EnvCorrsAM,'BinEdges',linspace(-0.5,1,30),'FaceColor','k')
histogram(EnvCorrs,linspace(-0.5,1,30),'FaceColor','g')
axis square
title('AM (black) vs Speech (green)')

suptitle('Correlations between stimuli')

% set(gcf,'PaperPosition',
orient(gcf,'landscape')
print(gcf,fullfile(savedir,'StimCorrs'),'-dpdf','-bestfit')


%% 












