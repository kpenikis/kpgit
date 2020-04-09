function FullTrialCTTS
% FullTrialCTTS
%
%   Random trial from each cell.
%   BIN: n spikes, n cells
% 
% KP, 2020-04
%
% close all
rng('shuffle')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
switch whichStim
    case 'AC'
        Stimuli      = 1:8;
        PickTrials   = 'all';
    case 'Speech'
        Stimuli      = 1:6;
        PickTrials   = 'sim';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minTr        = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binwidth     = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data files
fn = set_paths_directories;
savedir = fullfile(fn.figs,'PopDynamics');

switch whichStim
    case {'AC' 'DB'}
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
                
    case 'Speech'
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
end

q=load(fullfile(savedir,'Data',['CTTS_' whichStim '_' PickTrials]));
TD   = q.TD;
CTTS = q.Cell_Time_Trial_Stim;
ETTS = q.Env_Time_Trial_Stim;
clear q


%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%%
% Filter to identify cells with at least 8 trials
nTrialMat = nan(size(CTTS,1),size(CTTS,4));
for ist = 1:size(CTTS,4)
    CT  = permute(sum(CTTS(:,500+(1:TD(ist)),:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

% CellTypes
RS = UnitInfo.TroughPeak>0.43;
NS = UnitInfo.TroughPeak<=0.43;

% GCs = (sum(nTrialMat(:,Stimuli)<minTr,2)==0);
% theseCells = find(GCs & RS);



%% 1) Plot temporal relationship between peakRate and following envMax


hf_pReM=figure; 
set(hf_pReM,'Position',fullscreen)
hold on

for ist = Stimuli
    
    GCs = (sum(nTrialMat(:,ist)<minTr,2)==0);
    theseCells = find(GCs);
    
    % Find timestamps of all landmarks in envelope
    Env = mean(mean(ETTS(theseCells,1:(1000+TD(ist)),:,ist),3,'omitnan'),1);
    [allTS,~,~,~,minRatio] = find_peakRate(Env,1000);
    allTS = allTS(:,500+(1:TD(ist)));
    
    %         figure;
    %         plot(allTS(1,:),'k','LineWidth',1)
    %         hold on
    
    pRs = find(allTS(6,:))';
    eMs = nan(size(pRs));
    for ipr = 1:numel(pRs)
        try
            nextMax = find(allTS(4,(pRs(ipr)+1):min(pRs(ipr)+500,end)),1,'first');
        catch
            nextMax = find(allTS(4,pRs(ipr):end),1,'last');
        end
%         plot([pRs(ipr) eMs(ipr)],[allTS(4,eMs(ipr)) allTS(4,eMs(ipr))],'LineWidth',2)
        if ~isempty(nextMax)
            eMs(ipr) = pRs(ipr) + nextMax;
        elseif isempty(nextMax) && ipr<numel(pRs)
            keyboard
        end
    end
    
    subplot(2,numel(Stimuli)/2,ist)
    histogram(diff([pRs eMs],1,2),0:10:500,'FaceColor','k')
%     ylim([0 20])
end


% JUST sentences, concatenated
Stim = [3 4];
% Concatenate envelopes
Env = [];
for ist = Stim
    GCs = (sum(nTrialMat(:,ist)<minTr,2)==0);
    theseCells = find(GCs);
    Env = [Env mean(mean(ETTS(theseCells,500+(1:TD(ist)),:,ist),3,'omitnan'),1)];
end

[allTS,~,~,~,~] = find_peakRate(Env,1000);

pRs = find(allTS(6,:))';
eMs = nan(size(pRs));
for ipr = 1:numel(pRs)
    try
        nextMax = find(allTS(4,pRs(ipr)+(1:200)),1,'first');
    catch
        nextMax = find(allTS(4,pRs(ipr):end),1,'last');
    end
    %         plot([pRs(ipr) eMs(ipr)],[allTS(4,eMs(ipr)) allTS(4,eMs(ipr))],'LineWidth',2)
    if ~isempty(nextMax)
        eMs(ipr) = pRs(ipr) + nextMax;
    elseif isempty(nextMax) && ipr<numel(pRs)
        keyboard
    end
end

figure;
ih(1)=histogram(diff([pRs eMs],1,2),0:10:150,'FaceColor','k');
hold on
ih(2)=histogram(diff(pRs),0:10:500,'FaceColor','b');
xlim([0 500])
legend(ih,{'Time from peakRate until maxEnv' 'Time between peakRate events'})
title([whichStim ' ' num2str(Stim)])

print_eps_kp(gcf,fullfile(savedir,sprintf('%s_HistEventTiming',whichStim)))



%allTS = [ Env; ...
%     diff_Env;...
%     minEnv;...
%     peakEnv;...
%     minRate;...
%     peakRate];




%% 2) Pick one random trial from each cell, simulating single trial population

keyboard

nt = 100;
for ist = Stimuli
    
    GCs = (sum(nTrialMat(:,ist)<minTr,2)==0);
    theseCells = find(GCs & RS);

    Nbins = floor(TD(ist)/binwidth);
    BINS  =  0:binwidth:(binwidth*Nbins);
    XVALS =  BINS(1:end-1) + binwidth/2;
    
    
    nSpkBins  = nan(nt,Nbins);
    nCellsBin = nan(nt,Nbins);
    
    for it = 1:nt
        
        thisTr = ceil(rand(size(nTrialMat(:,ist))).*nTrialMat(:,ist));
        
        TrEnvs    = nan(numel(theseCells),TD(ist));
        raster    = nan(numel(theseCells),TD(ist));
%         peakRates = nan(numel(theseCells),TD(ist));
        for iu = 1:numel(theseCells)
            raster(iu,:)    = CTTS(theseCells(iu),500+(1:TD(ist)),thisTr(theseCells(iu)),ist);
            TrEnvs(iu,:)    = ETTS(theseCells(iu),500+(1:TD(ist)),thisTr(theseCells(iu)),ist);
%             [allTS,~, peakRates(iu,:),~,minRatio] = find_peakRate(TrEnvs(iu,:), 1000);
        end
        
        rasterbin = raster(:,1:BINS(end));
        nSpkBins(it,:) = sum(reshape(sum(rasterbin,1),[binwidth Nbins]),1);
        nCB = zeros(size(rasterbin,1),Nbins);
        for iu = 1:size(rasterbin,1)
            nCB(iu,sum(reshape(rasterbin(iu,:),[binwidth Nbins]),1)>0) = 1;
        end
        nCellsBin(it,:) = sum(nCB,1);
        
    end %it
    
    figure;
    plot(15+5.*(mean(mean(ETTS(:,500+(1:BINS(end)),:,ist),1,'omitnan'),3,'omitnan')./max(mean(mean(ETTS(:,500+(1:BINS(end)),:,ist),1,'omitnan'),3,'omitnan'))),...
        'Color',0.6*[1 1 1],'LineWidth',2)
%     plot(sum(raster,1),'k')
    hold on
    plot(XVALS,mean(nSpkBins,1))
    plot(XVALS,mean(nCellsBin,1))
    
    figure;
    plot((mean(mean(ETTS(:,500+(1:BINS(end)),:,ist),1,'omitnan'),3,'omitnan')./max(mean(mean(ETTS(:,500+(1:BINS(end)),:,ist),1,'omitnan'),3,'omitnan'))),...
        'Color',0.6*[1 1 1],'LineWidth',2)
%     plot(1+10.*mean(mean(ETTS(:,500+(1:BINS(end)),:,ist),1,'omitnan'),3,'omitnan'))
    hold on
    plot(XVALS,mean(nSpkBins,1)./mean(nCellsBin,1),'k','LineWidth',2)
    
end %ist



%% 3) For each cell, correlate PSTH to stimulus 

axmax = 1;
plim  = 0.01;

hf_scat=figure; 
set(hf_scat,'Position',fullscreen)
hold on

for ist = Stimuli
    
    GCs = (sum(nTrialMat(:,ist)<minTr,2)==0);
    theseCells = find(GCs & RS);
    
    r_a = nan(numel(theseCells),1);
    p_a = nan(numel(theseCells),1);
    r_dp = nan(numel(theseCells),1);
    p_dp = nan(numel(theseCells),1);
    r_dn = nan(numel(theseCells),1);
    p_dn = nan(numel(theseCells),1);
    
    for iu = 1:numel(theseCells)
        
        PSTH   = conv(mean(CTTS(theseCells(iu),:,:,ist),3,'omitnan'),convwin);
        PSTH   = PSTH(520+(1:TD(ist)));
        STIM   = mean(ETTS(theseCells(iu),500+(1:TD(ist)),:,ist),3,'omitnan');
        
        [r_a(iu),p_a(iu)] = corr(PSTH',STIM');
        [r_dp(iu),p_dp(iu)] = corr(PSTH(2:end)',diff(STIM)');
        
%         dp = diff(STIM');
%         dp(dp<0) = 0;
%         dn = -diff(STIM');
%         dn(dn<0) = 0;
%         
%         [r_dp(iu),p_dp(iu)] = corr(PSTH(2:end)',dp);
%         [r_dn(iu),p_dn(iu)] = corr(PSTH(2:end)',dn);
    end
    
    % Correct nans
    r_a(isnan(r_a)) = 0;
    p_a(isnan(p_a)) = 1;
    r_dp(isnan(r_dp)) = 0;
    p_dp(isnan(p_dp)) = 1;    
    
    
    AmpOnlyPos = (p_a<plim  & p_dp>=plim & r_a>0);
    AmpOnlyNeg = (p_a<plim  & p_dp>=plim & r_a<0);
    DrvOnlyPos = (p_a>=plim & p_dp<plim  & r_dp>0);
    DrvOnlyNeg = (p_a>=plim & p_dp<plim  & r_dp<0);
    Both_Q1    = (p_a<plim  & p_dp<plim & r_a>0 & r_dp>0);
    Both_Q2    = (p_a<plim  & p_dp<plim & r_a<0 & r_dp>0);
    Both_Q3    = (p_a<plim  & p_dp<plim & r_a<0 & r_dp<0);
    Both_Q4    = (p_a<plim  & p_dp<plim & r_a>0 & r_dp<0);
    
%     AmpOnly = (p_a<plim  & p_d>=plim);
%     DrvOnly = (p_a>=plim & p_d<plim);
%     Both    = (p_a<plim  & p_d<plim);
    
    figure(hf_scat); 
    hold on
    subplot(2,numel(Stimuli)/2,find(ist==Stimuli))
    
    plot([0 0],axmax*[-1 1],'-k')
    hold on
    plot(axmax*[-1 1],[0 0],'-k')
    plot(r_a,r_dp,'ok')
    plot(r_a(AmpOnlyPos|AmpOnlyNeg),r_dp(AmpOnlyPos|AmpOnlyNeg),'o','MarkerFaceColor','m','MarkerEdgeColor','k')
    plot(r_a(DrvOnlyPos|DrvOnlyNeg),r_dp(DrvOnlyPos|DrvOnlyNeg),'o','MarkerFaceColor','g','MarkerEdgeColor','k')
    plot(r_a(Both_Q1|Both_Q2|Both_Q3|Both_Q4),r_dp(Both_Q1|Both_Q2|Both_Q3|Both_Q4),'o','MarkerFaceColor','b','MarkerEdgeColor','k')
    plot(median(r_a),-1,'^','MarkerFaceColor','k','MarkerEdgeColor','none')
    plot(-1,median(r_dp),'>','MarkerFaceColor','k','MarkerEdgeColor','none')
    axis square
    xlabel('corr to amplitude')
    ylabel('corr to derivative')
    xlim(axmax*[-1 1])
    ylim(axmax*[-1 1])
    [r,p] = corr(r_a,r_dp);
    title(sprintf('%s %i\nr=%0.3f, p=%0.3f',whichStim,ist,r,p))
    
%     print_eps_kp(gcf,fullfile(savedir,sprintf('%s_CorrToStim_RS_%i',whichStim,ist)))
    
    
%     figure;
%     set(gcf,'Position',halfscreen./2)
%     subplot(1,2,1)
%     histogram(r_a,-axmax:0.05:axmax,'FaceColor','m')
%     hold on
%     plot(median(r_a),0,'^','MarkerFaceColor','k','MarkerEdgeColor','none')
%     xlim([-axmax axmax])
%     ylim([0 40])
%     title('Correlation to Amplitude')
%     
%     subplot(1,2,2)
%     histogram(r_dp,-axmax:0.05:axmax,'FaceColor','g')
%     hold on
%     plot(median(r_dp),0,'^','MarkerFaceColor','k','MarkerEdgeColor','none')
%     xlim([-axmax axmax])
%     ylim([0 40])
%     title('Correlation to Derivative')
    
%     print_eps_kp(gcf,fullfile(savedir,sprintf('%s_CorrHists_RS_%i',whichStim,ist)))
    
    
    MeanEnv = mean(mean(ETTS(theseCells,500+(1:TD(ist)),:,ist),3,'omitnan'),1);
%     [r,p] = corr(MeanEnv(2:end)',abs(diff(MeanEnv)'))
    
    
    % mean PSTHs of each type
% %     figure;
% %     plot(mean(mean(ETTS(theseCells,500+(1:TD(ist)),:,ist),3,'omitnan'),1),'k','LineWidth',2)
% %     hold on
% %     
% %     % Amplitude
% %     plot(conv(mean(mean(CTTS(theseCells(AmpOnlyPos),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(AmpOnlyPos),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         '-m','LineWidth',1)
% %     plot(conv(mean(mean(CTTS(theseCells(AmpOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(AmpOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         '-k','LineWidth',1)
% %     plot(conv(mean(mean(CTTS(theseCells(AmpOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(AmpOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         ':m','LineWidth',1)
% %     
% %     % Derivative
% % %     plot(conv(mean(mean(CTTS(theseCells(DrvOnlyPos),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin),...
% % %         '-r','LineWidth',1)
% % %     plot(conv(mean(mean(CTTS(theseCells(DrvOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin),...
% % %         '-k','LineWidth',1)
% % %     plot(conv(mean(mean(CTTS(theseCells(DrvOnlyNeg),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin),...
% % %         ':r','LineWidth',1)
% %     
% %     % Both
% %     plot(conv(mean(mean(CTTS(theseCells(Both_Q1),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(Both_Q1),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         'c','LineWidth',1)
% %     plot(conv(mean(mean(CTTS(theseCells(Both_Q2),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(Both_Q2),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         'b','LineWidth',1)
% %     plot(conv(mean(mean(CTTS(theseCells(Both_Q3),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(Both_Q3),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         'g','LineWidth',1)   
% %     plot(conv(mean(mean(CTTS(theseCells(Both_Q4),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)./max(conv(mean(mean(CTTS(theseCells(Both_Q4),500+(1:TD(ist)),:,ist),3,'omitnan'),1),convwin)),...
% %         'k','LineWidth',1) 
% %     
% %     xlim([0 1937])
    
end %ist


print_eps_kp(hf_scat,fullfile(savedir,sprintf('%s_CorrToStim_RS_all_shift20',whichStim)))




end



