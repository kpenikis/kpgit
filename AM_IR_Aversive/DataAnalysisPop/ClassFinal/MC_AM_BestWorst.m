function MC_AM_BestWorst
% MasterClass (all parameters defined at top of file)
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality.
%
%
%  KP, 2019-12
%

close all

varPar       = 'Full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'BestWorst';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
% Dur          = [10:20:150 200 300 400 500];
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials  = 'rand';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 500;
KernelType   = 'linear';
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
nBest        = 50;

rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'ClassAM','RawData');

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'CTTS_AM')); %Cell_Time_Trial_Stim
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% Load SU classification results
q = load(fullfile(fn.figs,'ClassAM',whichStim,varPar,'each','CR_each.mat'));
CReach = q.CR;
clear q


%%
% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
smallscreen = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];


%% Prepare to parse data

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


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ii = 1:numel(WinEnds)
    
    AnWin = WinBeg(ii):WinEnds(ii);
    fprintf('Dur: %i ms\n',WinEnds(ii)-WinBeg(ii)+1)
    
    %     for iTrSh = 1:numel(PickTrialss)
    %
    %         PickTrials = PickTrialss{iTrSh};
    
    % Define cells and stimuli
    [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
        'each', nTrialMat, UnitData,...
        theseStim, iRS, iNS, minTrs, convwin, AnWin );
    
    ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);
    
    
    % Find indices of best X (after "all" filtered)
    [dps,iSUdps] = sort(CReach.dprime,'descend');
    BestCells = iSUdps(1:nBest);
    
    % Count total N spikes for these N cells (average per trial)
    nSpk_Best = round(sum(sum(sum(mean(CTTS(BestCells,:,:,:),3,'omitnan'),2,'omitnan'),1)));
    
    % Now find X worst cells, that have roughly the same number of
    % spikes
    nSpk_Worst = 0;
    while nSpk_Worst<nSpk_Best
        if nSpk_Worst == 0
            iSU = length(iSUdps);
        else
            iSU = iSU-1;
        end
        nSpk_Worst = round(sum(sum(sum(mean(CTTS(iSU,:,:,:),3,'omitnan'),2,'omitnan'),1))) + nSpk_Worst;
    end
    
    WorstCells = iSUdps(iSU:length(iSUdps));
    
    
    
    % Train and test classifier
    
    % Best cells
    [S_AssMat_Best,E_AssMat_Best] = runSVMclass_SnE( CTTS(BestCells,:,:,:), ETTS(BestCells,:,:,:), ...
        BootstrapN, nStim, Dur, length(BestCells), PickTrials,TrainSize, TestSize, KernelType );
    
    % Worst cells
    [S_AssMat_Worst,E_AssMat_Worst] = runSVMclass_SnE( CTTS(WorstCells,:,:,:), ETTS(WorstCells,:,:,:), ...
        BootstrapN, nStim, Dur, length(WorstCells), PickTrials,TrainSize, TestSize, KernelType );
    
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %##########################################################################
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    keyboard
    
    %% Plot SU dprime vs N spikes
    
%     nSpks = (mean(sum(mean(CTTS,3,'omitnan'),2,'omitnan'),4));
%     
%     figure;
%     subplot(1,4,1)
%     hold on
%     
%     % Manually make boxplots
%     q5  = quantile(CReach.dprime,0.05);
%     q25 = quantile(CReach.dprime,0.25);
%     q75 = quantile(CReach.dprime,0.75);
%     q95 = quantile(CReach.dprime,0.95);
%     
%     plot([1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
%     fill(1+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
%     
%     ylabel('dprime')
%     set(gca,'Color','none')
%     box off
%     ylim([-0.5 2.5])
%     plotSpread(CReach.dprime,'distributionIdx',ones(size(CReach.dprime)),'distributionColors','k','showMM',3)
%     
%     
%     subplot(1,4,2:4)
%     plot(nSpks,CReach.dprime,'k.')
%     ylim([-0.5 2.5])
%     set(gca,'ytick',[],'Color','none')
%     box off
%     xlabel('N spks')
%     %         title('SU dprimes vs mean N spks per stim')
%     
%     savename = sprintf('SU_dps_%s_%s',varPar,whichStim);
%     print_eps_kp(gcf,fullfile(savedir,savename))
    
    
    %% Plot NEURAL results
    
    ConfMat = mean(S_AssMat_Best,3);
    muPC    = mean(diag(ConfMat))*100;
    dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
    
    ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);
    
    % Plot
    hf1(ii) = figure;
    imagesc(ConfMat)
    axis square
    caxis([-1 1])
    %         cmocean('ice') %ice
    cmocean('curl','pivot',0)
    colorbar
    ylabel('True stim')
    xlabel('Assigned')
    set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
    
    title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s (N=%i)',...
        muPC,dprime,KernelType,whichStim,whichCells,nUns))
    
    
    % Save figure
    savedir = fullfile(fn.figs,'ClassAM',whichStim,varPar);
    if ~exist(fullfile(savedir,'backupTables'),'dir')
        mkdir(fullfile(savedir,'backupTables'))
    end
    
    %         savename = sprintf('Res_v%s_Train%i_%s_%s_%s',varPar,...
    %             TrainSize,whichIrr,whichCells,PickTrials);
    savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s%i',varPar,WinEnds(ii)-WinBeg(ii)+1,TrainSize,tau,whichStim,'Best');
%     savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s',varPar,WinEnds(ii)-WinBeg(ii)+1,TrainSize,tau,whichStim,whichCells);
    
    %         print(hf(ii),fullfile(savedir,savename),'-dpdf')
    
    
    %% Save results to master table
    
    mastertablesavename = sprintf('CR_v%s_%s',varPar,whichCells);
    %         mastertablesavename = sprintf('CR_v%s_%s',varPar,whichIrr);
    thistablesavename   = savename;
    
    % Calculate performance for Env data
    ConfMat  = mean(E_AssMat,3);
    muPC_E   = mean(diag(ConfMat))*100;
    dprime_E = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
    
    
    CR1 = table;
    CR1.figname  = {savename};
    CR1.Stim     = {whichStim};
    CR1.Cells    = {whichCells};
    CR1.iC       = 0;
    CR1.Results  = {S_AssMat};
    CR1.PC       = muPC;
    CR1.dprime   = dprime;
    CR1.WinBeg   = WinBeg(ii);
    CR1.WinEnd   = WinEnds(ii);
    CR1.trials   = {PickTrials};
    CR1.rcorrOpt = {'norm'};
    CR1.nTemp    = {PSTHsize};
    CR1.nTrain   = TrainSize;
    CR1.nTest    = TestSize;
    CR1.conv     = {'exp'};
    CR1.tau      = tau;
    CR1.BsN      = BootstrapN;
    CR1.ResEnv   = {E_AssMat};
    CR1.PC_E     = muPC_E;
    CR1.dprime_E = dprime_E;
    
    
    % Load saved table
    clear q;
    try
        q=load(fullfile(savedir,mastertablesavename));
    end
    
    if ii>1 && ~exist('q','var')
        keyboard
    end
    
    if ~exist('q','var')
        
        % Save new table
        CR = CR1;
        save(fullfile(savedir,mastertablesavename),'CR','-v7.3')
        save(fullfile(savedir,'backupTables',thistablesavename),'CR1','-v7.3')
        
    else % Save updated table
        
        % Concatenate new data
        CR = q.CR;
        CR = [CR; CR1];
        
        save(fullfile(savedir,mastertablesavename),'CR','-v7.3')
        save(fullfile(savedir,'backupTables',thistablesavename),'CR1','-v7.3')
    end
    
    
    %     end % iTrSh
end %vary classification parameter


keyboard


% AnDur
dpe = [CR.dprime_E];
dpe(isinf(dpe)) = 6;

figure;
set(gcf,'Position',smallscreen)

subplot(2,1,1);
plot([CR.WinEnd]-500,dpe,'-m','LineWidth',2)
hold on
plot([CR.WinEnd]-500,[CR.dprime],'-k','LineWidth',2)
xlabel('Time from onset (ms)')
ylabel('d''')

subplot(2,1,2);
plot([CR.WinEnd]-500,[CR.dprime]-dpe,'-k','LineWidth',2)
xlabel('Time from onset (ms)')
ylabel('d'' Neural-Env')
ylim([-5 0])

print_eps_kp(gcf,fullfile(savedir,mastertablesavename))




% Shuffled vs sim trials, vary AnDur
figure;
hold on
ishf = find(strcmp(CR.trials,'rand'));
isim = find(strcmp(CR.trials,'sim'));
[~,iCR]=sort([CR.WinEnd]);
ip(1)=plot(CR.WinEnd(iCR(ismember(iCR,ishf)))-500,[CR.dprime(iCR(ismember(iCR,ishf)))],'o-r');
ip(2)=plot(CR.WinEnd(iCR(ismember(iCR,isim)))-500,[CR.dprime(iCR(ismember(iCR,isim)))],'o-k');
xlabel('Analysis window end')
ylabel('d''')
legend({'shuffled' 'sim.'},'Location','southeast')
title(whichCells)
print_eps_kp(gcf,fullfile(savedir,mastertablesavename))


% Shuffled vs sim trials
figure;
plot([0 3],[0 3],'-k')
hold on
scatter([CR.dprime(2:2:end)],[CR.dprime(1:2:end)],10.*[26 25 17 17 17 13],'k','filled')
xlabel('d'' simultaneous trials')
ylabel('d'' shuffled trials')
axis square

print_eps_kp(gcf,fullfile(savedir,mastertablesavename))


keyboard

figure; hold on
[~,iCR]=sort([CR.WinEnd]);
[~,iCRFR]=sort([CRFR.WinEnd]);
ip(1)=plot(CR.WinEnd(iCR)-500,[CR.dprime(iCR)],'o-b');
ip(2)=plot([CRFR.WinEnd(iCRFR)]-500,[CRFR.dprime(iCRFR)],'o-k');
xlabel('Analysis window end (ms)')
ylabel('d''')
legend(ip,{'temporal' 'rate'})
title('Classification in increasing windows')

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClassRcorr');
print_eps_kp(gcf,fullfile(savedir,'CR_vAnWin_50_TemporalRate'))

end
