function MC_cellpairs
% MasterClass (all parameters defined at top of file)
%  Can process AM or Speech data.
% 
%  MUST RE-RUN SPEECH!! wrong units labeled RS first time around
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality.
%
%
%  KP, 2019-01
%

close all

varPar       = 'Pairs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'Mar28-AM'; % BestRS SkipBestRS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = {'rand' 'sim'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 500;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
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

% Check that matching data files were imported
if size(Cell_Time_Trial_Stim,1)~=numel(UnitData)
    keyboard
end
if size(Cell_Time_Trial_Stim,1)<size(CReach,1)
    keyboard
end


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
figsavedir = fullfile(rootdir,whichStim,varPar,whichCells);
if ~exist(fullfile(figsavedir,'backupTables'),'dir')
    mkdir(fullfile(figsavedir,'backupTables'))
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
    
    % Get UnitData indices for finding SU dps
    [~,SU_UnNs,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,...
    theseStim, iRS, iNS, minTrs, convwin, AnWin );
    
    % Define cells and stimuli
    [CTTS,theseCells,~,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
        whichCells, nTrialMat, UnitData,...
        theseStim, iRS, iNS, minTrs, convwin, AnWin );
    
    [aa,iCReach] = ismember(theseCells,SU_UnNs);
    if ~all(aa)
        keyboard
    end
    
    ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);
    
    for iC1 = 1:(numel(theseCells)-1)
        for iC2 = (iC1+1):numel(theseCells)
            
            UseCells = theseCells([iC1 iC2]);
            isRS     = UnitInfo(UseCells,:).TroughPeak>0.43;
            
            for TrM = 1:numel(PickTrials)
                
                TrMethod = PickTrials{TrM};
                
                % Train and test classifier
                [S_AssMat_NC,~,~] = runSVMclass_SnE( CTTS([iC1 iC2],:,:,:), ETTS([iC1 iC2],:,:,:), ...
                    BootstrapN, nStim, Dur, length(UseCells), TrMethod,TrainSize, TestSize, KernelType );
                
                %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                %##########################################################################
                %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                
                
                %% Plot NEURAL results
                
                ConfMat = mean(S_AssMat_NC,3);
                
                muPC    = mean(diag(ConfMat))*100;
                dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
                
                ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);
                
                % Plot
                hf(ii) = figure;
                imagesc(ConfMat)
                axis square
                caxis([-1 1])
                cmocean('curl','pivot',0)
                colorbar
                ylabel('True stim')
                xlabel('Assigned')
                set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
                
                title(sprintf('%0.1f%%, d''=%0.2f\nSVM (%s)  |  %s (%i & %i)',...
                    muPC,dprime,whichStim,whichCells,UseCells(1),UseCells(2)))
                
                
                % Save figure
%                 savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s_%i',varPar,WinEnds(ii)-WinBeg(ii)+1,TrainSize,tau,whichStim,whichCells,nCell);
                savename = sprintf('Res_%s-%i_Train%i_conv%i_%s_%s-%i-%i_%s',varPar,WinEnds(ii)-WinBeg(ii)+1,TrainSize,tau,whichStim,whichCells,UseCells(1),UseCells(2),TrMethod);
                
                print(hf(ii),fullfile(figsavedir,savename),'-dpdf')
                
                
                %% Save results to master table
                
                mastertablesavename = sprintf('CR_v%s_%s',varPar,whichCells);
                thistablesavename   = savename;
                
                CR1 = table;
                CR1.figname  = {savename};
                CR1.Stim     = {whichStim};
                CR1.Cells    = {whichCells};
                CR1.iC       = 2;
                CR1.Results  = {S_AssMat_NC};
                CR1.PC       = muPC;
                CR1.dprime   = dprime;
                CR1.WinBeg   = WinBeg(ii);
                CR1.WinEnd   = WinEnds(ii);
                CR1.trials   = {TrMethod};
                CR1.rcorrOpt = {'norm'};
                CR1.nTemp    = {PSTHsize};
                CR1.nTrain   = TrainSize;
                CR1.nTest    = TestSize;
                CR1.conv     = {'exp'};
                CR1.tau      = tau;
                CR1.iUnData  = {UseCells};
                CR1.iCReach  = {iCReach([iC1 iC2])};
                CR1.isRS     = {isRS};
                
                
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
                
            end % TrM
            close all
        end% iC2
    end % iC1
end %vary classification parameter

iCpairs = [];
for iC1 = 1:(numel(theseCells)-1)
    for iC2 = (iC1+1):numel(theseCells)
        iCpairs = [iCpairs; [iC1 iC2] ] ;
    end
end

pcdiffs = CR.PC(2:2:end)-CR.PC(1:2:end);
dpdiffs = CR.dprime(2:2:end)-CR.dprime(1:2:end);

figure;
subplot(2,1,1)
histogram(dpdiffs)
subplot(2,1,2)
histogram(pcdiffs)

% Better simultaneous: NS & RS cell, opposite response timing
figure; plot(mean(Cell_Time_Trial_Stim(110,:,:,3),3,'omitnan')); hold on; plot(mean(Cell_Time_Trial_Stim(121,:,:,3),3,'omitnan'));
% Better shuffled: RS & RS cell, similar response timing
figure; plot(mean(Cell_Time_Trial_Stim(102,:,:,3),3,'omitnan')); hold on; plot(mean(Cell_Time_Trial_Stim(121,:,:,3),3,'omitnan'));


% Does principle hold for top 5 cells of each?  Not really.. 
[dpds,idpds] = sort(dpdiffs);
ShuffBest = idpds(1:5);
SimBest   = idpds((end-4):end);

corr(mean(CTTS(iCpairs(ShuffBest(1),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(ShuffBest(1),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(ShuffBest(2),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(ShuffBest(2),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(ShuffBest(3),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(ShuffBest(3),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(ShuffBest(4),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(ShuffBest(4),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(ShuffBest(5),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(ShuffBest(5),2),:,:,7),3,'omitnan')')

corr(mean(CTTS(iCpairs(SimBest(1),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(SimBest(1),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(SimBest(2),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(SimBest(2),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(SimBest(3),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(SimBest(3),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(SimBest(4),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(SimBest(4),2),:,:,7),3,'omitnan')')
corr(mean(CTTS(iCpairs(SimBest(5),1),:,:,7),3,'omitnan')',mean(CTTS(iCpairs(SimBest(5),2),:,:,7),3,'omitnan')')


keyboard

% Plot ranked SU vs subpops
pcr_SUvPools

% Plot min max PC and d' as a function of N cells
pcr_minmaxPC

% Plot PC and d' for each stimulus as a function of N cells
pcr_byStim

% Plot SU dprime vs N spikes
pcr_SUdpFR




% AnDur
dpe = [CR.dprime_E];
dpe(isinf(dpe)) = 6;

figure;
set(gcf,'Position',tallsmall)

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



end
