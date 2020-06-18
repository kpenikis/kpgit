function PopResp_CTTS
% 
% PopResp_CTTS
% 
% 2020-03, updated 2020-04 

close all

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
% SORT
sortBy       = 'peakFRquant_Lat'; % 'FRrange'; 
if strcmp(whichStim,'AC')
    sortStim     = 3;
elseif strcmp(whichStim,'Speech')
    sortStim     = 4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHING
% convtype     = 'gauss';
% tau          = 20;
% window       = gausswin(tau);
% window       = window-min(window);
% convwin      = window/sum(window);

convtype     = 'exp';   
tau          = 10;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
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
q = load(fullfile(rootdir,whichStim,'Full','each','CR_each.mat'));
CReach = q.CR;
clear q


%%
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
widehalf    = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];
widenarrow  = [1 scrsz(4)/4 scrsz(3) scrsz(4)/4];

% Set figsavedir
figsavedir = fullfile(fn.figs,'PopResp','NSmin');
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
% [CTTS,theseUns,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
%     'PopResp', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, AnWin, convtype );
% 
% if size(CTTS,1) ~= size(Cell_Time_Trial_Stim,1)
%     keyboard
% end
% FR_vec = permute(mean(CTTS,3,'omitnan'),[1 2 4 3]);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
minTrs = 12;

% Get UnitData indices for matching to CR each
[CTTS,theseUns,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, AnWin, 'exp' );

FR_vec = permute(mean(CTTS,3,'omitnan'),[1 2 4 3]).*1000;

% Load encoding time estimation
% load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/avgPropUp.mat')

% New RS/NS labels
flagRS = find(UnitInfo(theseUns,:).TroughPeak>0.43);
flagNS = find(UnitInfo(theseUns,:).TroughPeak<=0.43);

% Get envelope data
ETTS = Env_Time_Trial_Stim(theseUns,AnWin,:,theseStim);


%%
% rangeFR = mean([max(FR_vec(:,:,2),[],2)-min(FR_vec(:,:,2),[],2) ...
%     max(FR_vec(:,:,3),[],2)-min(FR_vec(:,:,3),[],2) ...
%     max(FR_vec(:,:,4),[],2)-min(FR_vec(:,:,4),[],2) ],2);
% figure;
% plot(rangeFR(theseUns),CReach.dprime,'ok')
% [r,p]=corr(rangeFR(theseUns),CReach.dprime)
% 
% 
% peakFR = mean([max(FR_vec(:,:,2),[],2) ...
%     max(FR_vec(:,:,3),[],2) ...
%     max(FR_vec(:,:,4),[],2) ],2);
% figure;
% plot(peakFR(theseUns),CReach.dprime,'ok')
% [r,p]=corr(peakFR(theseUns),CReach.dprime)
% 
% figure;
% plot(peakFR,rangeFR,'ok')
% [r,p]=corr(peakFR,rangeFR)


%% Sort data

clear i_sorted
%

switch sortBy
    case 'FRrange'
        
        ytickset = 0;
        yticklab = {''};
        
        % RS cells
        rangeFR    = mean([max(FR_vec(flagRS,:,2),[],2)-min(FR_vec(flagRS,:,2),[],2) ...
            max(FR_vec(flagRS,:,3),[],2)-min(FR_vec(flagRS,:,3),[],2) ...
            max(FR_vec(flagRS,:,4),[],2)-min(FR_vec(flagRS,:,4),[],2) ],2);
        Qbounds    = quantile(rangeFR,[0 0.2 0.4 0.6 0.8 1]);
        [~,ipkRS]  = max(FR_vec(flagRS,:,sortStim),[],2);
        
        figure; hold on
        
        dataRS  = [];
        isortRS = [];
        for iq=1:5
            idx = find(rangeFR>=Qbounds(iq) & rangeFR<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkRS(idx),'descend');
            isortRS  = [isortRS; idx(idx_sort)];
            dataRS   = [dataRS; rangeFR(idx(idx_sort)) ipkRS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab ['RS' num2str(iq)]];
            
            subplot(2,3,iq)
            plot(avgPropUp(intersect(theseUns,flagRS(idx)) ), ...
                CReach.dprime( (ismember(theseUns,flagRS(idx))) ),'ok')
            title(['Q' num2str(iq)])
            ylim([-0.2 3])
            xlim([0 0.7])
            axis square
            [r,p] = corr(avgPropUp(intersect(theseUns,flagRS(idx)) ), ...
                CReach.dprime( (ismember(theseUns,flagRS(idx))) ))
        end
        
        % NS cells
        rangeFR    = mean([max(FR_vec(flagNS,:,2),[],2)-min(FR_vec(flagNS,:,2),[],2) ...
            max(FR_vec(flagNS,:,3),[],2)-min(FR_vec(flagNS,:,3),[],2) ...
            max(FR_vec(flagNS,:,4),[],2)-min(FR_vec(flagNS,:,4),[],2) ],2);
        Qbounds    = quantile(rangeFR,[0 1]);
        [~,ipkNS]  = max(FR_vec(flagNS,:,sortStim),[],2);
        
        dataNS  = [];
        isortNS = [];
        for iq=1
            idx = find(rangeFR>=Qbounds(iq) & rangeFR<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkNS(idx),'descend');
            isortNS  = [isortNS; idx(idx_sort)];
            dataNS   = [dataNS; rangeFR(idx(idx_sort)) ipkNS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab 'NS'];
        end
        
        sortdata     = [dataRS; dataNS];
        i_sorted     = [flagRS(isortRS); flagNS(isortNS)];
        ytickset     = cumsum(ytickset);
        
    
    case 'peakFRquant_Lat'
        
        ytickset = 0;
        yticklab = {''};
        
        %/\/\/\/\/\/\/\/\
        %    RS cells
        %\/\/\/\/\/\/\/\/
%         [pkRS,ipkRS] = max(ThisData(flagRS,:,sortStim),[],2);
        [pkRS,isort]  = rankPeakFR(FR_vec(flagRS,:,:));        % to group
        [~,isback] = sort(isort);
        pkRS = pkRS(isback);
        if strcmp(whichStim,'AC') && sortStim==3                 % to sort
            [~,ipkRS] = max(FR_vec(flagRS,1:250,sortStim),[],2);
        elseif strcmp(whichStim,'Speech') 
            [~,ipkRS] = max(FR_vec(flagRS,:,sortStim),[],2);
        else
            keyboard
            [~,ipkRS] = max(FR_vec(flagRS,:,sortStim),[],2);
        end
        Qbounds      = quantile(pkRS,[0 0.2 0.4 0.6 0.8 1]);
        
        lats_RS = nan(ceil(numel(flagRS)/5),5);
        dataRS  = [];
        isortRS = [];
        for iq=1:5
            idx = find(pkRS>=Qbounds(iq) & pkRS<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkRS(idx),'descend');
            isortRS  = [isortRS; idx(idx_sort)];
            dataRS   = [dataRS; pkRS(idx(idx_sort)) ipkRS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab ['RS' num2str(iq)]];
            
            lats_RS(1:size(lats,1),iq) = lats';
        end
        
        %/\/\/\/\/\/\/\/\
        %    NS cells
        %\/\/\/\/\/\/\/\/
%         [pkNS,ipkNS] = max(ThisData(flagNS,:,sortStim),[],2);
        [pkNS,isort]  = rankPeakFR(FR_vec(flagNS,:,:));        % to group
        [~,isback] = sort(isort);
        pkNS = pkNS(isback);
        if strcmp(whichStim,'AC') && sortStim==3                 % to sort
%             [~,ipkNS] = max(FR_vec(flagNS,1:250,sortStim),[],2);
            [~,ipkNS] = min(FR_vec(flagNS,1:250,sortStim),[],2);
        elseif strcmp(whichStim,'Speech') 
            [~,ipkNS] = min(FR_vec(flagNS,:,sortStim),[],2);
        else
            keyboard
            [~,ipkNS] = min(FR_vec(flagNS,:,sortStim),[],2);
        end
        Qbounds      = quantile(pkNS,[0 1]);
        
        lats_NS = nan(numel(flagNS),1);
        dataNS  = [];
        isortNS = [];
        for iq=1
            idx = find(pkNS>=Qbounds(iq) & pkNS<=Qbounds(iq+1));
            [lats,idx_sort] = sort(ipkNS(idx),'descend');
            isortNS  = [isortNS; idx(idx_sort)];
            dataNS   = [dataNS; pkNS(idx(idx_sort)) ipkNS(idx(idx_sort))];
            ytickset = [ytickset size(lats,1)];
            yticklab = [yticklab 'NS'];
            
            lats_NS(1:size(lats,1),1) = lats';
        end
        
        sortdata     = [dataRS; dataNS];
        i_sorted     = [flagRS(isortRS); flagNS(isortNS)];
        ytickset     = cumsum(ytickset);
end

% % 
% % %% Plot latency distributions
% % 
% % % % % set(gcf,'Position',tallhalf)
% % % % if strcmp(whichStim,'AC')
% % % %     set(gcf,'Position',fullscreen)
% % % % elseif strcmp(whichStim,'Speech')
% % % %     set(gcf,'Position',widehalf)
% % % % else
% % % %     keyboard
% % % % end
% % % % 
% % % % 
% % % % subplot(1,size(ThisData,3),sortStim)
% % 
% % 
% % hf2 = figure;       
% % set(gca,'ColorOrder',cmocean('thermal',6))
% % hold on
% % plot(nan,nan,'.')
% % for iq = 1:size(lats_RS,2)
% %     histogram(lats_RS(:,iq),0:500,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2)
% % end
% % histogram(lats_NS,0:500,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k');
% % 
% % ylim([0 1])
% % xlim([0 500])
% % set(gca,'Color','none')
% % box off
% % grid on
% % 
% % % Stats
% % 
% % % 
% % % keyboard
% % % 
% % % savestr = sprintf('%s_%s_st%i_%s%i',whichStim,sortBy,sortStim,convtype,tau);
% % % print_eps_kp(hf2,fullfile(figsavedir,['LatDistr_' savestr]))
% % 
% % 
% % %% Plot results

hf1 = figure;
% set(gcf,'Position',tallhalf)
if strcmp(whichStim,'AC')
    set(gcf,'Position',fullscreen)
elseif strcmp(whichStim,'Speech')
    set(gcf,'Position',widehalf)
else
    keyboard
end

ThisData = FR_vec;

for ist = 1:size(FR_vec,3)
    
    subplot(1,size(FR_vec,3),ist)
    thisRespPlot
    
end


%%
savestr = sprintf('%s_%s_st%i_%s%i',whichStim,sortBy,sortStim,convtype,tau);
% print_eps_kp(hf1,fullfile(figsavedir,['PopResp_' savestr]))

% set(hf1,'PaperOrientation','landscape')
% print(hf1,'-dpdf','-r500','-fillpage', fullfile(figsavedir,['PopResp_' whichStim '_' sortBy '_v2']))

% set(hf1,'PaperOrientation','portrait')
% print(hf1,'-dpdf','-r500','-fillpage', fullfile(figsavedir,['PopResp_' whichStim '_' sortBy '_portrait']))



%% Also: 
% grand avg PSTH (RS cells) 
% correlate to amplitude and derivative



% first, convert amplitude to dB! 
% then, MinPeakProminence can be +6 (doubling of loudness) 



% close all
hf2 = figure;
% set(gcf,'Position',widenarrow)
set(gcf,'Position',widehalf)
for ist = 1:size(FR_vec,3)
    
    % Get variables
    Env = mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1);
    Env = Env./max(Env);
    
    LdB = sound2dB_AM(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1),75); %mode(TrialData.SPL)
    
    Drv = diff(Env);
    Drv(Drv<0) = 0;
    Drv = Drv./max(Drv);
    GAH = mean(FR_vec(flagRS,:,ist),1)./max(mean(FR_vec(flagRS,:,ist),1));
    
    % Yulia's method
    [allTS, ~, peakRate, peakEnv, minRatio] = find_peakRate(LdB, 1000);
    
    
    % Plot data
    subplot(3,size(FR_vec,3),ist)
%     plot(Drv,'r','LineWidth',2)
    hold on
%     plot(Env,'b','LineWidth',2)
    plot(allTS(1,:),'c','LineWidth',2)
    plot(GAH,'k','LineWidth',2)
    plot(find(peakRate),GAH(find(peakRate)),'.r','MarkerSize',20)
    plot(find(peakEnv),GAH(find(peakEnv)),'.g','MarkerSize',20)
    plot(find(allTS(3,:)),GAH(find(allTS(3,:))),'.c','MarkerSize',20)
    ylim([0 1])
%     ylim([0 0.25])
    title(ist)
    
    % Calculate correlations
    [C_env,lags] = xcorr(GAH,Env,60);
    [C_drv,lags] = xcorr(GAH(2:end),Drv,60);
    
    subplot(3,size(FR_vec,3),size(FR_vec,3)+ist)
    plot(lags,C_drv,'r','LineWidth',2)
    hold on
    plot(lags,C_env,'b','LineWidth',2)
    
%     % Yulia's method
%     subplot(3,size(FR_vec,3),size(FR_vec,3)*2+ist)
%     
end




%% Population sparseness

keyboard


hf2 = figure;
set(gcf,'Position',fullscreen)

for ist = 1:size(FR_vec,3)
    
    subplot(5,size(FR_vec,3),ist)
%     plot([0 500],[0.6135 0.6135],'k')
    hold on
    
    PopSps = nan(1,size(FR_vec,2));
    PopSps_NS = nan(1,size(FR_vec,2));
    for ims = 1:size(FR_vec,2)
        PopSps(1,ims) = calculateSparseness(FR_vec(flagRS,ims,ist));
        PopSps_NS(1,ims) = calculateSparseness(FR_vec(flagNS,ims,ist));
    end
    plot(PopSps_NS,'b','LineWidth',2)
    plot(PopSps,'m','LineWidth',2)
    ylim([0 0.9])
    set(gca,'Color','none')
    box off
    
    subplot(5,size(FR_vec,3),size(FR_vec,3)+ist)
    plot(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1,'omitnan'),'k','LineWidth',2)
    if strcmp(whichStim,'AC')
        ylim([0 0.9])
    else
        ylim([0 0.02])
    end
    set(gca,'Color','none','xtick',[],'ytick',[])
    box off
end

print_eps_kp(hf2,fullfile(figsavedir,['PopSps_' savestr]))


end
