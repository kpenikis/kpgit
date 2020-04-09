function LatDur_CTTS
% 
% PopResp_CTTS
% 

global dpSU AnWin whichLat FRdiffThresh

close all


whichClass   = 'ActVec';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'top10'; %'SpecRest'; %'Q_pkFR'; %'dpRank_RS'; % pkFR_RS
whichLat     = 'pk_time'; 'onset'; 
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
% Convolution
tau          = 10;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crAmt        = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution for latency histogram (and envelope?)
tau          = 20;
%gaussian window
convwin2     = gausswin(tau);
convwin2     = convwin2-min(convwin2);
convwin2     = convwin2/sum(convwin2);
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
ETTS  = permute(mean(mean(q.Env_Time_Trial_Stim,3,'omitnan'),1,'omitnan'),[4 2 1 3]);
clear q

% Load SU classification results
q = load(fullfile(rootdir,whichStim,whichClass,'each','CR_each.mat'));
CReach = q.CR;
clear q

% % Load Quantile classification results
% q=load(fullfile(rootdir,whichStim,whichClass,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
% CR_Qpfr_Sp = q.CR;
% clear q

% Load encoding time estimation
% load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/avgPropUp.mat')


%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',12)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallhalf    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
widehalf    = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

% Set figsavedir
figsavedir = fullfile(fn.figs,'PopResp','LatDur');
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
minTrs = 12;

% Get UnitData indices for matching to CR each
[CTTS,theseUns,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, 1:size(Cell_Time_Trial_Stim,2), 'exp' );

ETTS = ETTS(theseStim,AnWin);

if size(CTTS,1)~=size(CReach,1)
    keyboard
end

% New RS/NS labels
flagRS = find(UnitInfo(theseUns,:).TroughPeak>0.43);
flagNS = find(UnitInfo(theseUns,:).TroughPeak<=0.43);

% Flag non-significant SU dps
% UnSig = bootstrap4significance(CReach);

nStim = size(CReach(1,:).Results{:},1);

FRdiffThresh = 0.1;
% pdms = 500;
DurMax = 350;
upd_lat = makedist('uniform','Lower',1,'Upper',500);
upd_dur = makedist('uniform','Lower',1,'Upper',DurMax);


if strcmp(whichCells,'top10')
    
    hfLD=figure;
    set(hfLD,'Position',widehalf)
    hold on
    
    hfLP=figure;
    set(hfLP,'Position',widehalf)
    hold on
    
    hfLatHist=figure;
    set(hfLatHist,'Position',widehalf./[1 1 1 2])
    
    hfPSTH=figure;
    set(hfPSTH,'Position',widehalf./[1 1 1 2])
    
    
    % Get all SU stim d'
    dpSU = nan(size(CReach,1),nStim);
    for ii=1:size(CReach,1)
        dpSU(ii,:) = dp_from_ConfMat(mean(CReach(ii,:).Results{:},3,'omitnan'),crAmt);
    end
    
    % [ Sp/Rest lat dur d' ]
    Result = nan(size(CReach,1),5,8);
    ResultStats = struct;
    
    LatHistRS_Cat = [];
    LatHistNS_Cat = [];
    PSTH_RS       = [];
    PSTH_NS       = [];
    
    r_AmpRS = nan(8,1);    p_AmpRS = nan(8,1);
    r_DrvRS = nan(8,1);    p_DrvRS = nan(8,1);
    r_AmpNS = nan(8,1);    p_AmpNS = nan(8,1);
    r_DrvNS = nan(8,1);    p_DrvNS = nan(8,1);
    
    for ist = 1:8
        
%         iC_Spec = find(dpSU(flagRS,ist)>prctile(dpSU(flagRS,ist),whatPerc));
%         iC_Rest = find(dpSU(flagRS,ist)<=prctile(dpSU(flagRS,ist),whatPerc));
        
        [foo,ifoo] = sort(dpSU(flagRS,ist),'descend');
        iC_Spec = ifoo(1:10);
        iC_Rest = ifoo(11:end);
        
        %------------------------------------------------------------------
        %%----Rest (RS non-specialist cells)
        %------------------------------------------------------------------
        Result = getLatDur(flagRS(iC_Rest),ist,CTTS,Result,CReach);
        
        % Lat vs dur (size peak height)
        figure(hfLD); hold on
        subplot(2,4,ist);
        hold on
        scatter(Result(flagRS(iC_Rest),2,ist), Result(flagRS(iC_Rest),3,ist), 2*Result(flagRS(iC_Rest),4,ist),...
            'o','MarkerFaceColor',[1 1 1]*0.04,'MarkerEdgeColor','none')
        
        % Lat vs peak FR (size d')
        figure(hfLP); hold on
        subplot(2,4,ist);
        hold on
        scatter(Result(flagRS(iC_Rest),2,ist), Result(flagRS(iC_Rest),4,ist), ceil(40*dpSU(flagRS(iC_Rest),ist))+20,...
            'o','MarkerFaceColor',[1 1 1]*0.04,'MarkerEdgeColor','none')        
        
        
        %------------------------------------------------------------------
        %%%----Specialists
        %------------------------------------------------------------------
        Result = getLatDur(flagRS(iC_Spec),ist,CTTS,Result,CReach);
        
        % Lat vs dur (size peak height)
        figure(hfLD); hold on
        subplot(2,4,ist);
        hold on
        scatter(Result(flagRS(iC_Spec),2,ist), Result(flagRS(iC_Spec),3,ist), 2*Result(flagRS(iC_Spec),4,ist),...
            'o','MarkerFaceColor','m','MarkerEdgeColor','none')
        
        % Lat vs peak FR (size d')
        figure(hfLP); hold on
        subplot(2,4,ist);
        hold on
        scatter(Result(flagRS(iC_Spec),2,ist), Result(flagRS(iC_Spec),4,ist), ceil(40*dpSU(flagRS(iC_Spec),ist))+20,...
            'o','MarkerFaceColor','m','MarkerEdgeColor','none')
        
        
        %------------------------------------------------------------------
        %%%---- NS
        %------------------------------------------------------------------
        Result = getLatDur(flagNS,ist,CTTS,Result,CReach);
        
        
        %------------------------------------------------------------------
        % Finish plots
        %------------------------------------------------------------------
        figure(hfLD); hold on
        subplot(2,4,ist);
        hold on
        xlim([0 500])
        ylim([0 DurMax])
        set(gca,'Color','none')
        switch whichLat
            case 'pk_time'
                xlabel('Time of peak FR (ms)')
            case 'onset'
                xlabel('Time exceeds 1/4 of max FR (ms)')
        end
        ylabel('Duration above 1/4 of max FR (ms)')
        
        figure(hfLP); hold on
        subplot(2,4,ist);
        hold on
        xlim([0 500])
        ylim([0 150])
        set(gca,'Color','none')
        switch whichLat
            case 'pk_time'
                xlabel('Time of peak FR (ms)')
            case 'onset'
                xlabel('Time exceeds 1/4 of max FR (ms)')
        end
        ylabel('Peak FR (sp/s)')
        
        
        %------------------------------------------------------------------
        % Histogram of peak times
        %------------------------------------------------------------------
        
        figure(hfLatHist); hold on
        
        subplot(2,4,ist);
        hold on
        histogram(Result(flagRS(iC_Rest),2,ist),1:10:500,'FaceColor',[1 1 1]*0.03,'EdgeColor','none','FaceAlpha',1)
        histogram(Result(flagNS,2,ist),1:10:500,'FaceColor','g','EdgeColor','none','FaceAlpha',1)
        histogram(Result(flagRS(iC_Spec),2,ist),1:10:500,'FaceColor','m','EdgeColor','none','FaceAlpha',1)
        set(gca,'Color','none')
        xlim([0 500])
        ylim([0 20])
        switch whichLat
            case 'pk_time'
                xlabel('Time of peak FR (ms)')
            case 'onset'
                xlabel('Time exceeds 1/4 of max FR (ms)')
        end
        
        
        %------------------------------------------------------------------
        % Grand avg PSTHs
        %------------------------------------------------------------------
        
        figure(hfPSTH); hold on
        
        subplot(2,4,ist);
        hold on
        plot(mean(mean(CTTS(flagRS,AnWin,:,ist),3,'omitnan'),1)*1000,...
            'Color',0.08*[1 1 1],'LineWidth',2)
        plot(mean(mean(CTTS(flagNS,AnWin,:,ist),3,'omitnan'),1)*1000,...
            'g','LineWidth',2)
        set(gca,'Color','none')
        xlim([0 500])
        ylim([0 20])
        
        
        PSTH_RS = [PSTH_RS mean(mean(CTTS(flagRS,AnWin,:,ist),3,'omitnan'),1)];
        PSTH_NS = [PSTH_NS mean(mean(CTTS(flagNS,AnWin,:,ist),3,'omitnan'),1)];
        
        
        %------------------------------------------------------------------
        % Statistics
        %------------------------------------------------------------------
        
        [~,p_lat_RS] = kstest(Result(flagRS,2,ist),'cdf',upd_lat);
        [~,p_lat_NS] = kstest(Result(flagNS,2,ist),'cdf',upd_lat);
        
        ResultStats.StimLat_Uniform_RS(ist) = p_lat_RS;
        ResultStats.StimLat_Uniform_NS(ist) = p_lat_NS;
        
        fprintf('Stim %i: uniform latency dist\n  RS p=%0.1e\n  NS p=%0.1e\n',...
            ist,p_lat_RS,p_lat_NS)
        
        
        %------------------------------------------------------------------
        % Relate to amplitude
        %------------------------------------------------------------------
%         find_peakRate(ETTS(ist,:), 1000)
        
        % RS cells
        
        LatHistRS = histcounts(Result(flagRS,2,ist),0.5:500.5);
        LatHistRS = conv(LatHistRS,convwin2,'same');
        LatHistRS_Cat = [LatHistRS_Cat LatHistRS];
        
        [r_AmpRS(ist),p_AmpRS(ist)] = corr(LatHistRS',ETTS(ist,:)');
        [r_DrvRS(ist),p_DrvRS(ist)] = corr(LatHistRS(1:end-1)',diff(ETTS(ist,:))');
        
        % NS cells
        LatHistNS = histcounts(Result(flagNS,2,ist),0.5:500.5);
        LatHistNS = conv(LatHistNS,convwin2,'same');
        LatHistNS_Cat = [LatHistNS_Cat LatHistNS];
        
        [r_AmpNS(ist),p_AmpNS(ist)] = corr(LatHistNS',ETTS(ist,:)');
        [r_DrvNS(ist),p_DrvNS(ist)] = corr(LatHistNS(1:end-1)',diff(ETTS(ist,:))');
        
    end %ist
    
    
    %% Correlate peak times to envelope
    
    Env_Cat = reshape(ETTS',[numel(ETTS) 1]);
    
    figure; hold on
    plot(LatHistRS_Cat./max(LatHistRS_Cat),'k','LineWidth',2)
    plot(LatHistNS_Cat./max(LatHistNS_Cat),'g','LineWidth',2)
    plot(Env_Cat./max(Env_Cat),':c','LineWidth',2)
    plot(diff(Env_Cat)./max(diff(Env_Cat)),':b','LineWidth',2)
    
    [r_AmpRS_Lat,p_AmpRS_Lat] = corr(LatHistRS_Cat',Env_Cat)
    [r_DrvRS_Lat,p_DrvRS_Lat] = corr(LatHistRS_Cat(1:end-1)',diff(Env_Cat))
    [r_Dv2RS_Lat,p_Dv2RS_Lat] = corr(LatHistRS_Cat(1:end-2)',diff(Env_Cat,2));
    
    [r_AmpNS_Lat,p_AmpNS_Lat] = corr(LatHistNS_Cat',Env_Cat)
    [r_DrvNS_Lat,p_DrvNS_Lat] = corr(LatHistNS_Cat(1:end-1)',diff(Env_Cat))
    [r_Dv2NS_Lat,p_Dv2NS_Lat] = corr(LatHistNS_Cat(1:end-2)',diff(Env_Cat,2));
    
    % Now use avg PSTHs
    [r_AmpRS_PSTH,p_AmpRS_PSTH] = corr(PSTH_RS',Env_Cat)
    [r_DrvRS_PSTH,p_DrvRS_PSTH] = corr(PSTH_RS(1:end-1)',diff(Env_Cat))
    [r_Dv2RS_PSTH,p_Dv2RS_PSTH] = corr(PSTH_RS(1:end-2)',diff(Env_Cat,2));
    
    [r_AmpNS_PSTH,p_AmpNS_PSTH] = corr(PSTH_NS',Env_Cat)
    [r_DrvNS_PSTH,p_DrvNS_PSTH] = corr(PSTH_NS(1:end-1)',diff(Env_Cat))
    [r_Dv2NS_PSTH,p_Dv2NS_PSTH] = corr(PSTH_NS(1:end-2)',diff(Env_Cat,2));
    
    
    
    %% Save figures
    
    savename = sprintf('%s_%s_%s_LatDur',whichStim,whichCells,whichLat);
    print_eps_kp(hfLD,fullfile(figsavedir,savename))
    
    savename = sprintf('%s_%s_%s_LatPeak',whichStim,whichCells,whichLat);
    print_eps_kp(hfLP,fullfile(figsavedir,savename))
    
    savename = sprintf('%s_%s_%s_LatHIST',whichStim,whichCells,whichLat);
    print_eps_kp(hfLatHist,fullfile(figsavedir,savename))
    
    savename = sprintf('%s_PSTHs',whichStim);
    print_eps_kp(hfPSTH,fullfile(figsavedir,savename))
    
end

keyboard

% Find sound events
%     allTS = [TS; ...
%              diff_loudness;...
%              minEnv;...
%              peakEnv;...
%              minRate;...
%              peakRate];
allTS = find_peakRate(stim_AM, 1000, [1 length(stim_AM)]./1000, 'broadband');



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% X change TopEach to Top10 cells
% X add histogram of latencies
% X add NS cells)
% X check stats

% X grand PSTH RS vs NS

%   *relate to stimulus amplitude
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



end





