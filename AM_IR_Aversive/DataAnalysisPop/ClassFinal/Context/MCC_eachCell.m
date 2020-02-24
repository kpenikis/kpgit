function MCC_eachCell
% MasterClass (all parameters defined at top of file)
%
%  CONTEXT CLASSIFICATION
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality. 
%
%
%  KP, 2020-01
%

close all

varPar       = 'OnlyRate';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'each'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;%[10:20:150 200 300 400 500];
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'rand';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 300;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
winlen       = 500;
convwin      = exp(-lambda*(1:winlen));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        rootdir = fullfile(fn.figs,'ClassContext','AM');
        rawdata = 'CTTSC_AM_sim.mat';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassContext','Speech');
        rawdata = 'CTTSC_Speech_sim.mat';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
        k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
        kfns = fieldnames(k);
        SegDurs  = structfun(@length,k);
        
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
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];

% Set figsavedir
figsavedir = fullfile(rootdir,varPar,whichCells);
if ~exist(fullfile(figsavedir,'backupTables'),'dir')
    mkdir(fullfile(figsavedir,'backupTables'))
end


%% Prepare to parse data

theseStim  = 1:size(Cell_Time_Trial_Stim,5);

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for iSeg = 1:size(Cell_Time_Trial_Stim,4)
    
    AnWin = 1:SegDurs(iSeg);
    fprintf('%s\n',kfns{iSeg})
    
    nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,5));
    for ist = 1:size(Cell_Time_Trial_Stim,5)
        CT  = permute(sum(Cell_Time_Trial_Stim(:,AnWin,:,iSeg,ist),2),[1 3 5 2 4]);
        nTrialMat(:,ist) = sum(~isnan(CT),2);
    end
    
        
    % For rate only classifier: shuffle spiketimes within trial
    if strcmp(varPar,'OnlyRate')
        CTTSC = shuffleSpikeTimes(Cell_Time_Trial_Stim(:,AnWin,:,iSeg,:));
    else
        CTTSC = Cell_Time_Trial_Stim(:,AnWin,:,iSeg,:);
    end
    
    
    % Define cells and stimuli
    [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( CTTSC, ...
        whichCells, nTrialMat, UnitData,theseStim, iRS, iNS, minTrs, convwin, 1:size(CTTSC,2) );
        
    % Step through each unit
    for ii = 1:numel(theseCells)
        
        iUn = theseCells(ii);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %##########################################################################
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        DataIn = single(CTTS(ii,:,:,:));
        
        % Train and test classifier
        S_AssMat = runSVMclass_SnE( DataIn, DataIn, ...
            BootstrapN, nStim, Dur, size(DataIn,1), PickTrials, TrainSize, TestSize, KernelType );
        
%         S_AssMat = runSVMclassFR( DataIn, BootstrapN, nStim, Dur, size(DataIn,1), PickTrials, TrainSize, TestSize, KernelType );
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %##########################################################################
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        %% Plot NEURAL result
        
        ConfMat = mean(S_AssMat,3);
        muPC    = mean(diag(ConfMat))*100;
        dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
        
        ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);
        
        % Plot
        hf(ii) = figure;
        imagesc(ConfMat)
        axis square
        caxis([-1 1])
        %         cmocean('ice') %ice
        cmocean('curl','pivot',0)
        colorbar
        ylabel('True Context')
        xlabel('Assigned')
        set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
        
        title(sprintf('%0.1f%%, d''=%0.2f\n%s  |  %s (iUn=%i)',...
            muPC,dprime,kfns{iSeg},whichCells,iUn))
        
        
        % Save figure
        
        savename = sprintf('ResContext_%s_%s-%i_Train%i_conv%i_%s',kfns{iSeg},whichCells,iUn,TrainSize,tau,whichStim);
        
        print(hf(ii),fullfile(figsavedir,savename),'-dpdf')
        
        
        %% Save results to master table
        
        mastertablesavename = sprintf('CR_%s',whichCells);
        thistablesavename   = savename;
        
        % Load data into table
        CR1 = table;
        CR1.figname  = {savename};
        CR1.Stim     = {whichStim};
        CR1.Cells    = {whichCells};
        CR1.iC       = iUn;
        CR1.Seg      = iSeg;
        CR1.SegStr   = kfns(iSeg);
        CR1.Results  = {S_AssMat};
        CR1.PC       = muPC;
        CR1.dprime   = dprime;
        CR1.trials   = {PickTrials};
        CR1.rcorrOpt = {'norm'};
        CR1.nTemp    = {PSTHsize};
        CR1.nTrain   = TrainSize;
        CR1.nTest    = TestSize;
        CR1.conv     = {'exp'};
        CR1.tau      = tau;
        CR1.BsN      = BootstrapN;
        
        
        % Load saved table
        clear q;
        if ~exist('CR','var')
        try
            q=load(fullfile(figsavedir,mastertablesavename));
            CR = q.CR;
        end
        end
        
        if ~exist('CR','var')
            
            % Save new table
            CR = CR1;
            save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
            save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
            
        else % Save updated table
            
            % Concatenate new data
            CR = [CR; CR1];
            save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
            save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
        end
        
        if mod(iUn,10)==0
            close all
        end
    end %iUn
end %vary classification parameter


keyboard



figure;
plot(CReach.dprime,CR.dprime,'k.')
hold on
plot([0 3],[0 3])
axis square

signrank(CReach.dprime,CR.dprime)
median(CReach.dprime-CR.dprime)


iUns = unique(CReach.iC);

dpHat_temp = nan(numel(iUns),1);
for iiu = 1:numel(iUns)
    icr = CR.iC==iUns(iiu);
    dpHat_temp(iiu) = mean(CR.dprime(icr));
end


iUns = unique(CReach.iC);

dpHat_full = nan(numel(iUns),1);
for iiu = 1:numel(iUns)
    icr = CReach.iC==iUns(iiu);
    dpHat_full(iiu) = mean(CReach.dprime(icr));
end


dpContext = sort(dpHat_full,'descend');
figure;
plot(1:numel(dpContext),dpContext,'.k')
hold on

CReach = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/Speech/Full/each/CR_each.mat');
CReach = CReach.CR;

dpShape = sort(CReach.dprime,'descend');

plot(1:numel(dpShape),dpShape,'.b')



keyboard


pcr_SUdpFR


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

print_eps_kp(gcf,fullfile(figsavedir,mastertablesavename))




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
print_eps_kp(gcf,fullfile(figsavedir,mastertablesavename))


% Shuffled vs sim trials
figure;
plot([0 3],[0 3],'-k')
hold on
scatter([CR.dprime(2:2:end)],[CR.dprime(1:2:end)],10.*[26 25 17 17 17 13],'k','filled')
xlabel('d'' simultaneous trials')
ylabel('d'' shuffled trials')
axis square

print_eps_kp(gcf,fullfile(figsavedir,mastertablesavename))


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
figsavedir = fullfile(fn.figs,'StimClassRcorr');
print_eps_kp(gcf,fullfile(figsavedir,'CR_vAnWin_50_TemporalRate'))

end
