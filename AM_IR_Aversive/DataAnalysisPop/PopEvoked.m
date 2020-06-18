function PopEvoked
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
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHING
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
figsavedir = fullfile(fn.figs,'PopResp','EventEvoked');
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
    'each', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, 1:size(Cell_Time_Trial_Stim,2), 'exp' );

FR_vec = permute(mean(CTTS,3,'omitnan'),[1 2 4 3]).*1000;

% Load encoding time estimation
% load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/avgPropUp.mat')

% New RS/NS labels
flagRS = find(UnitInfo(theseUns,:).TroughPeak>0.43);
flagNS = find(UnitInfo(theseUns,:).TroughPeak<=0.43);

% Get envelope data
ETTS = Env_Time_Trial_Stim(theseUns,:,:,theseStim);


%% Also: 
% grand avg PSTH (RS cells) 
% correlate to amplitude and derivative

PPT = 200;

AllSegs_Ons = [];
AllSegs_PkD = [];
AllSegs_PkE = [];
Tdiff_EtoD =[]; 
Tdiff_DtoE =[]; 
Tdiff_DtoD =[]; 
Tdiff_MtoD =[]; 

hf2 = figure;
set(gcf,'Position',widehalf)
hold on

hfst = figure;
set(gcf,'Position',widehalf)

hTD=figure;
hold on

AMrates = [2 4 8 16 32]; 

for ist = 2:6 % 1:size(FR_vec,3)  %[2 3 4 5 7 8]
    
    % First, convert amplitude to dB scale! 
    % then, MinPeakProminence can be +6 (a doubling of loudness)
    switch whichStim
        case 'AC'
            LdB = sound2dB_AM(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1),60); %mode(TrialData.SPL)
%             LdB = mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1);
        case 'Speech'
            LdB = sound2dB_Sp(mean(mean(ETTS(:,:,:,ist),3,'omitnan'),1)); %mode(TrialData.SPL)
    end
    
    % Yulia's method, adapted
    [allTS, ~, peakRate, peakEnv, minRatio] = find_peakRate(LdB, 1000);
    
    
    % Find events just near AnWin
    Ons_all = find(allTS(3,:));
    PkD_all = find(allTS(6,:));
    PkE_all = find(allTS(4,:));
    
    % Remove events outside of AnWin
    Ons = Ons_all;
    Ons(Ons<(AnWin(1)-10)) = [];
    Ons(Ons>AnWin(end))    = [];
    PkD = PkD_all;
    PkD(PkD<(AnWin(1)-10)) = [];
    PkD(PkD>AnWin(end))    = [];
    PkE = PkE_all;
    PkE(PkE<(AnWin(1)-10)) = [];
    PkE(PkE>AnWin(end))    = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Timing relationships
    %%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:numel(PkE)
        % Find PkD preceding and following
        diffs = PkD_all-PkE(ii);
        % All events
        Tdiff_EtoD = [Tdiff_EtoD diffs(abs(diffs)<=PPT)];
        % Just preceding and following
%         Dpr = max(diffs(diffs<0));
%         Dfl = min(diffs(diffs>0));
%         Tdiff_EtoD = [Tdiff_EtoD Dpr Dfl];
%         Tdiff_EtoD = [Tdiff_EtoD diffs(diffs==max(diffs(diffs<0))) diffs(diffs==min(diffs(diffs>0)))];
    end
    for ii = 1:numel(PkD)
        
        % Find PkE preceding and following
        diffs = PkE_all-PkD(ii);
        % All events
        Tdiff_DtoE = [Tdiff_DtoE diffs(abs(diffs)<=PPT)];
        % Just preceding and following
%         Epr = max(diffs(diffs<0));
%         Efl = min(diffs(diffs>0));
%         Tdiff_DtoE = [Tdiff_DtoE Epr Efl];
        
        % Find PkD following
        diffs = PkD_all-PkD(ii);
        Dfl = min(diffs(diffs>0));
        Tdiff_DtoD = [Tdiff_DtoD Dfl];
        
%         figure(hTD); hold on
%         plot(Efl,Dfl,'ok')
        
        
        % Find MinEnv preceding
        diffs = Ons_all-PkD(ii);
        % All events
        %         Tdiff_DtoE = [Tdiff_DtoE diffs(abs(diffs)<=PPT)];
        % Just preceding and following
        MEpr = max(diffs(diffs<0));
        %if AM
%         Tdiff_MtoD = [Tdiff_MtoD MEpr/(1000/AMrates(ist-1)) ];
        %if speech
        Tdiff_MtoD = [Tdiff_MtoD MEpr ];
        
    end
    
    % Create and save plot of envelope + landmarks
%     figure;
%     plot((AnWin(1)-10):AnWin(end),allTS(1,(AnWin(1)-10):AnWin(end)),'k','LineWidth',8)
% %     plot((AnWin(1)-10):AnWin(end),LdB(1,(AnWin(1)-10):AnWin(end)),'k','LineWidth',8)
%     hold on
%     try
%     plot([PkD; PkD],[-120*ones(size(PkD)); 0*ones(size(PkD))],':r','LineWidth',2);
%     ip(2)=plot(PkD,allTS(1,PkD),'.r','MarkerSize',50);
%     end
%     try
%     plot([PkE; PkE],[-120*ones(size(PkE)); 0*ones(size(PkE))],':g','LineWidth',2);
%     ip(3)=plot(PkE,allTS(1,PkE),'.g','MarkerSize',50);
%     end
%     try
%     plot([Ons; Ons],[-120*ones(size(Ons)); 0*ones(size(Ons))],':b','LineWidth',2);
%     ip(1)=plot(Ons,allTS(1,Ons),'.b','MarkerSize',50);
%     end
%     
%     title(ist)
%     xlim([AnWin(1)-10 AnWin(end)])
% %     ylim([-120 0])
%     ylim([0 1])
%     set(gca,'Color','none')
%     if ist==8
%         legend(ip,{'Onset' 'Peak Deriv.' 'Peak Env'},'Location','best')
%     end
    
%     print_eps_kp(gcf,fullfile(figsavedir,['Stim_' whichStim '_amp_st' num2str(ist)]))
    
    
    
    % Clip segment of FR around each event
    
    for ii = 1:numel(PkD)
        FRseg = [];
        FRseg = mean(FR_vec(flagRS,(-PPT:PPT)+PkD(ii),ist),1);
%         plot(-PPT:PPT,FRseg,'r')
        AllSegs_PkD = [AllSegs_PkD; FRseg];
    end
    for ii = 1:numel(PkE)
        FRseg = [];
        FRseg = mean(FR_vec(flagRS,(-PPT:PPT)+PkE(ii),ist),1);
%         plot(-PPT:PPT,FRseg,'g')
        AllSegs_PkE = [AllSegs_PkE; FRseg];
    end
    for ii = 1:numel(Ons)
        FRseg=[];
        FRseg = mean(FR_vec(flagRS,(-PPT:PPT)+Ons(ii),ist),1);
%         plot(-PPT:PPT,FRseg,'c')
        AllSegs_Ons = [AllSegs_Ons; FRseg];
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subplot for each stimulus
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    range(AllSegs_Ons((end-numel(Ons)+1):end,:)')
    
    figure(hfst);
    hold on
    subplot(2,4,ist)
    hold on
%     try
%     plot(-PPT:PPT,mean(AllSegs_PkD((end-numel(PkD)+1):end,:),1),'r','LineWidth',2)
%     end
    try
    plot(-PPT:PPT,mean(AllSegs_PkE((end-numel(PkE)+1):end,:),1),'g','LineWidth',2)
    end
    try
    plot(-PPT:PPT,mean(AllSegs_Ons((end-numel(Ons)+1):end,:),1),'b','LineWidth',2)
    end
    ylim([0 20])
    set(gca,'Color','none')
    
    
%     % Normalize signals before calculating correlations
%     Env = allTS(1,AnWin);
%     Env = Env-min(Env);
%     Env = Env./max(Env);
%     
%     Drv = allTS(2,AnWin);
% %     Drv = Drv-min(Drv);
% %     Drv(Drv<0) = 0;
%     Drv = abs(Drv);
%     Drv = Drv./max(Drv);
%     
%     GAH = mean(FR_vec(flagRS,AnWin,ist),1);
%     GAH = GAH-min(GAH);
%     GAH = GAH./max(GAH);
%     
%     %     figure;
%     %     plot(Env,'g','LineWidth',2)
%     %     hold on
%     %     plot(Drv,'r','LineWidth',2)
%     %     plot(GAH,'k','LineWidth',2)
%     
%     % Calculate correlations
%     [C_env,lags] = xcorr(GAH,Env,60,'coeff');
%     [C_drv,lags] = xcorr(GAH,Drv,60,'coeff');
%     
%     figure(hf2);
%     
%     subplot(3,size(FR_vec,3),ist)
%     plot(GAH,'k','LineWidth',2)
%     hold on
%     plot(Env,'Color',[0 0.7 0.1],'LineWidth',2)
%     plot(Drv,'Color',[0.7 0 0],'LineWidth',2)
%     xlim([0 500])
%     
%     subplot(3,size(FR_vec,3),size(FR_vec,3)+ist)
%     plot([0 0],[0 1],'k')
%     hold on
%     plot(lags,C_drv,'Color',[0.7 0 0],'LineWidth',2)
%     plot(lags,C_env,'Color',[0 0.7 0.1],'LineWidth',2)
%     [md,imd] = max(C_drv);
%     [me,ime] = max(C_env);
%     plot(lags(imd),md,'.k','MarkerSize',20)
%     plot(lags(ime),me,'.k','MarkerSize',20)
%     ylim([0 1])
end

% Save evoked each stimulus figure
% print_eps_kp(hfst,fullfile(figsavedir,['EvokedStim_' whichStim '-400']))
% 
% % Save correlations figure
% print_eps_kp(hf2,fullfile(figsavedir,['Correlations_' whichStim]))




% Plot extracted FR segments surrounding each type of event 

mu_PkD  = mean(AllSegs_PkD,1);
std_PkD = std(AllSegs_PkD,1);
mu_PkE  = mean(AllSegs_PkE,1);
std_PkE = std(AllSegs_PkE,1);
mu_Ons  = mean(AllSegs_Ons,1);
std_Ons = std(AllSegs_Ons,1);


hfE=figure;
hold on

% Plot Evoked FR std or confidence intervals
subplot(6,1,1:5)
hold on

plot(repmat(-PPT:PPT,[size(AllSegs_PkE,1) 1])',AllSegs_PkE','k','LineWidth',0.5)

%STD
% plot(-PPT:PPT,mu_PkD+std_PkD,'r:','LineWidth',2)
% plot(-PPT:PPT,mu_PkD-std_PkD,'r:','LineWidth',2)
plot(-PPT:PPT,mu_PkE+std_PkE,'g:','LineWidth',2)
plot(-PPT:PPT,mu_PkE-std_PkE,'g:','LineWidth',2)
% plot(-PPT:PPT,mu_Ons+std_Ons,'b:','LineWidth',2)
% plot(-PPT:PPT,mu_Ons-std_Ons,'b:','LineWidth',2)

% Add mean Evoked FR
% plot(-PPT:PPT,mu_PkD,'Color',[0.5 0 0],'LineWidth',4)
plot(-PPT:PPT,mu_PkE,'Color',[0 0.4 0.1],'LineWidth',4)
% plot(-PPT:PPT,mu_Ons,'Color',[0 0 0.5],'LineWidth',4)
xlim([-50 100])
ylim([4 18])
set(gca,'Color','none')
title(whichStim)


subplot(6,1,6)
plot([Tdiff_EtoD; Tdiff_EtoD],[0;1]+zeros(2,length(Tdiff_EtoD)),'k','LineWidth',1)
% histogram(Tdiff_EtoD,-PPT:10:PPT,'FaceColor','r')
% xlim([-PPT PPT])
xlim([-50 100])
set(gca,'Color','none')



hfD=figure;
hold on

% Plot Evoked FR std or confidence intervals
subplot(6,1,1:5)
hold on

plot(repmat(-PPT:PPT,[size(AllSegs_PkD,1) 1])',AllSegs_PkD','k','LineWidth',0.5)

%STD
plot(-PPT:PPT,mu_PkD+std_PkD,'r:','LineWidth',2)
plot(-PPT:PPT,mu_PkD-std_PkD,'r:','LineWidth',2)
% plot(-PPT:PPT,mu_PkE+std_PkE,'g:','LineWidth',2)
% plot(-PPT:PPT,mu_PkE-std_PkE,'g:','LineWidth',2)
% plot(-PPT:PPT,mu_Ons+std_Ons,'b:','LineWidth',2)
% plot(-PPT:PPT,mu_Ons-std_Ons,'b:','LineWidth',2)

% Add mean Evoked FR
plot(-PPT:PPT,mu_PkD,'Color',[0.5 0 0],'LineWidth',4)
% plot(-PPT:PPT,mu_PkE,'Color',[0 0.4 0.1],'LineWidth',4)
% plot(-PPT:PPT,mu_Ons,'Color',[0 0 0.5],'LineWidth',4)
xlim([-50 100])
ylim([4 18])
set(gca,'Color','none')
title(whichStim)

subplot(6,1,6)
plot([Tdiff_DtoE; Tdiff_DtoE],[0;1]+zeros(2,length(Tdiff_DtoE)),'k','LineWidth',1)
% histogram(Tdiff_DtoE,-PPT:2.5:PPT,'FaceColor','g')
% xlim([-PPT PPT])
xlim([-50 100])
set(gca,'Color','none')



hfONS=figure;
hold on

% Plot Evoked FR std or confidence intervals
subplot(6,1,1:5)
hold on

plot(repmat(-PPT:PPT,[size(AllSegs_Ons,1) 1])',AllSegs_Ons','k','LineWidth',0.5)

%STD
% plot(-PPT:PPT,mu_PkD+std_PkD,'r:','LineWidth',2)
% plot(-PPT:PPT,mu_PkD-std_PkD,'r:','LineWidth',2)
% plot(-PPT:PPT,mu_PkE+std_PkE,'g:','LineWidth',2)
% plot(-PPT:PPT,mu_PkE-std_PkE,'g:','LineWidth',2)
plot(-PPT:PPT,mu_Ons+std_Ons,'b:','LineWidth',2)
plot(-PPT:PPT,mu_Ons-std_Ons,'b:','LineWidth',2)

% Add mean Evoked FR
% plot(-PPT:PPT,mu_PkD,'Color',[0.5 0 0],'LineWidth',4)
% plot(-PPT:PPT,mu_PkE,'Color',[0 0.4 0.1],'LineWidth',4)
plot(-PPT:PPT,mu_Ons,'Color',[0 0 0.5],'LineWidth',4)
xlim([-50 100])
ylim([4 18])
set(gca,'Color','none')
title(whichStim)

subplot(6,1,6)
plot([-1*Tdiff_MtoD; -1*Tdiff_MtoD],[0;1]+zeros(2,length(Tdiff_MtoD)),'k','LineWidth',1)
% histogram(Tdiff_DtoE,-PPT:2.5:PPT,'FaceColor','g')
% xlim([-PPT PPT])
xlim([-50 100])
set(gca,'Color','none')


%%% SAVE FIGS
print_eps_kp(hfE,fullfile(figsavedir,['EvokedResp_' whichStim '-Env']))
print_eps_kp(hfD,fullfile(figsavedir,['EvokedResp_' whichStim '-Drv']))
print_eps_kp(hfONS,fullfile(figsavedir,['EvokedResp_' whichStim '-Ons']))



%%
keyboard
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
