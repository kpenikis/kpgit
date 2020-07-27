function MC_RandSubpop
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
%  KP, 2020-01; updated 2020-03 for ActVec; 2020-04 back to Projections
%

close all

whichClass   = 'Full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'Rand_NS'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'rand';
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
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TrainSize    = 11;
TestSize     = 1;
minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPools       = 500;
PoolSize     = 10;

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
figsavedir = fullfile(rootdir,whichStim,whichClass,whichCells,'Sum');
if ~exist(fullfile(figsavedir),'dir')
    mkdir(fullfile(figsavedir))
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
iNS = find(UnitInfo.TroughPeak<=0.43);

% For rate only classifier: shuffle spiketimes within trial
ShuffOpt = 0;
if strcmp(whichClass,'OnlyRate')
    ShuffOpt = 1;
    convwin  = 1;
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ii = 1:numel(WinEnds)
    
    AnWin = WinBeg(ii):WinEnds(ii);
    fprintf('Dur: %i ms\n',WinEnds(ii)-WinBeg(ii)+1)
    
    for fc = (1:nPools) %PoolStart
        
        FirstCell = fc;
        
        for nc = 1:numel(PoolSize)
            
            NumCells = PoolSize(nc);
            
            % Define cells and stimuli
            [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
                'each', nTrialMat, UnitData,...
                theseStim, iRS, iNS, minTrs, convwin, AnWin );
            
            ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);
            
            % Get indices of RS//NS cells (from just "theseCells")
            iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
            iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);
            
            
            % Set the subpopulation of cells to use
            switch whichCells
                
                case 'Rand_RS'
                    [pkFRsort,ipkFR] = rankPeakFR(CTTS(iRS,:,:,:));
                    iRands     = ipkFR(randperm(size(ipkFR,1),NumCells));
                    UseCells   = iRS(iRands);
                    SUdps      = CReach(UseCells,:).dprime
                    
                case 'Rand_NS'
                    [pkFRsort,ipkFR] = rankPeakFR(CTTS(iNS,:,:,:));
                    iRands     = ipkFR(randperm(size(ipkFR,1),NumCells));
                    UseCells   = iNS(iRands);
                    SUdps      = CReach(UseCells,:).dprime
            end
            
                        
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if strcmp(whichClass,'OnlyTemp')
                keyboard
            elseif strcmp(whichClass,'Nspk')
                S_AssMat = runSVMclassFR( CTTS(UseCells,:,:,:), BootstrapN, nStim, Dur, length(UseCells), PickTrials, TrainSize, TestSize, KernelType );
            elseif strcmp(whichClass,'Full')
                %Train and test classifier
                [S_AssMat,~,~] = runSVMclass_notNorm( sum(CTTS(UseCells,:,:,:),1), sum(ETTS(UseCells,:,:,:),1), ...
                    BootstrapN, nStim, Dur, 1, PickTrials,TrainSize, TestSize, KernelType, ShuffOpt );
            elseif strcmp(whichClass,'ActVec')
                S_AssMat = runSVM_ActVec(CTTS(UseCells,:,:,:),BootstrapN, nStim, Dur, length(UseCells), PickTrials, TrainSize, TestSize, KernelType, ShuffOpt );
            end
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %##########################################################################
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
            
            %% Plot NEURAL results
            
            ConfMat = mean(S_AssMat,3);
            
            muPC    = mean(diag(ConfMat))*100;
            dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
            
%             ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);
            
            % Plot
%             hf(ii) = figure;
%             imagesc(ConfMat)
%             axis square
%             caxis([-1 1])
%             cmocean('curl','pivot',0)
%             colorbar
%             ylabel('True stim')
%             xlabel('Assigned')
%             set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
%             
%             title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s (N=%i,%i)',...
%                 muPC,dprime,KernelType,whichStim,whichCells,NumCells,FirstCell))
            
            
            % Save figure
            savename = sprintf('Res_%s-%i_%s_%s_%i_%i',...
                whichClass,WinEnds(ii)-WinBeg(ii)+1,whichStim,whichCells,NumCells,FirstCell);
            
%             print(hf(ii),fullfile(figsavedir,savename),'-dpdf')
            
            
            %% Save results to master table
            
            mastertablesavename = sprintf('CR_%s_%s',whichClass,whichCells);
            thistablesavename   = savename;
            
            CR1 = table;
            CR1.figname  = {savename};
            CR1.Stim     = {whichStim};
            CR1.Cells    = {whichCells};
            CR1.iC       = FirstCell;
            CR1.nC       = NumCells;
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
            CR1.UnIDs    = {theseCells(UseCells)};
            CR1.CRids    = {UseCells};
            CR1.SUdps    = {SUdps};
%             CR1.TrRes    = {TrialResults};
            
            
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
%                 save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
                
            else % Save updated table
                
                % Concatenate new data
                CR = [CR; CR1];
                save(fullfile(figsavedir,mastertablesavename),'CR','-v7.3')
%                 save(fullfile(figsavedir,'backupTables',thistablesavename),'CR1','-v7.3')
            end
            
        end % nc
    end % fc
end %vary classification parameter


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
rootdir = fullfile(fn.figs,'StimClassRcorr');
print_eps_kp(gcf,fullfile(rootdir,'CR_vAnWin_50_TemporalRate'))

end
