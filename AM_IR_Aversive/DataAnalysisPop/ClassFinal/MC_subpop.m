function MC_subpop(minTrs,whichStim,whichClass)  %(exclSpec,exclNonSig)
% MasterClass (all parameters defined at top of file)
%  Can process AM or Speech data.
% 
%  Instead of feeding spiking data into SVM (or whatever classifier), input
%  just a scalar for each class comparison. This allows for independent
%  information from many neurons, without a drastic increase of
%  dimensionality.
%
%
%  KP, 2020-01
%

% close all

% whichClass   = 'Full';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'dpRank_RS'; %'Q_pkFR'; %'allRS'; %'dpRank_RS';  
exclNonSig   = 0;
exclSpec     = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
WinBeg       = 501 * ones(size(Dur));
% WinBeg       = [501 551 601 651 751 901];
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'rand';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
% whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 50;
KernelType   = 'linear';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSTHsize     = 'Train-1';
TestSize     = 1;
TrainSize    = minTrs - TestSize;
% minTrs       = TrainSize + TestSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(whichCells,'Q_pkFR')
    PoolStart    = [0.99999 0.8 0.6 0.4 0.2];
    PoolSize     = 0.2;
elseif strcmp(whichCells,'dpRank_RS') || strcmp(whichCells,'maxdp_RS')
    PoolStart    = 1; %
%     PoolSize     = [3 5 10 nan];%[1 5 10 20 50];
else
    PoolStart    = 1;
    PoolSize     = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng('shuffle')


%% Load data

fn = set_paths_directories;
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
% q = load(fullfile(rootdir,whichStim,whichClass,'each','CR_each.mat'));
try 
    q = load(fullfile(rootdir,whichStim,whichClass,['minTrs' num2str(minTrs)],'CR_each.mat'));
catch
    q = load(fullfile(rootdir,whichStim,'ActVec',['minTrs' num2str(minTrs)],'CR_each.mat'));
end
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
figsavedir = fullfile(rootdir,whichStim,whichClass,['minTrs' num2str(minTrs)],whichCells);
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
        theseStim  = 1:size(Cell_Time_Trial_Stim,4);
end

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);


% Flag non-significant SU dps
if exclNonSig
    UnSig = bootstrap4significance(CReach);
end

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
    
    % For rate only classifier: shuffle spiketimes within trial
%     if strcmp(whichClass,'OnlyRate')
%         Cell_Time_Trial_Stim(:,AnWin,:,:) = shuffleSpikeTimes(Cell_Time_Trial_Stim(:,AnWin,:,:));
%     end
    
    % Define cells and stimuli
    [CTTS,theseCells,~,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
        'each', nTrialMat, UnitData,theseStim, iRS, iNS, minTrs, convwin, AnWin );
    if numel(theseCells)~=size(CReach,1)
        keyboard
    end
    
    ETTS = Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);

    % Get indices of RS//NS cells (from just "theseCells")
    iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
    iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);
    
    %%%% Set PoolSize now, based on how many RS cells total with this minTrs
    if numel(iRS)>100
        PoolSize = [3 5 10 20 30 40 50 numel(iRS)];
    elseif numel(iRS)>50
        PoolSize = [3 5 7 10 15 20 30 numel(iRS)];
    elseif numel(iRS)>30
        PoolSize = [3 5 7 10 20 numel(iRS)];
    else
        PoolSize = [3 5 7 10 numel(iRS)];
    end
    
    for fc = 1:numel(PoolStart) %PoolStart
        
        if PoolStart(fc)<1
            FirstCell = round(numel(iRS)*PoolStart(fc));
        else
            FirstCell = PoolStart(fc);
        end
        
        for nc = 1:numel(PoolSize)
            
            if PoolSize(nc)<1
                NumCells = round(numel(iRS)*PoolSize(nc));
            else
                NumCells = PoolSize(nc);
            end
            
%             if (FirstCell+NumCells-1) > numel(iRS)
%                 fprintf(' SKIPPING ')
%                 continue
%             end
            
            % Set the subpopulation of cells to use
            switch whichCells
                
                case 'all'
                    UseCells   = 1:size(CTTS,1);
                    
                case 'allRS'
                    UseCells   = iRS;
                    
                case 'allNS'
                    UseCells   = iNS;
                    
                case 'pkFR_RS'
                    [pkFRsort,ipkFR] = rankPeakFR(CTTS(iRS,:,:,:));
                    UseCells   = iRS(ipkFR(FirstCell+(0:(NumCells-1))));
%                     SUdps      = CReach(UseCells,:).dprime
                    
                case 'dpRank_RS'
                    [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    if isnan(NumCells)
                        UseCells   = iRS(iSUdps);
                    else
                        UseCells   = iRS(iSUdps(FirstCell+(0:(NumCells-1))));
                    end
%                     SUdps      = CReach(UseCells,:).dprime
                
                case 'maxdp_RS'
                    [maxdps,iSUdps] = sort_maxdp(CReach(iRS,:));
                    UseCells   = iRS(iSUdps(FirstCell+(0:(NumCells-1))));
                    
                case 'Q_pkFR'
                    [pkFRsort,ipkFR] = rankPeakFR(CTTS(iRS,:,:,:));
                    UseCells   = iRS(ipkFR( FirstCell-(NumCells:-1:1)+1 ));
%                     SUdps      = CReach(UseCells,:).dprime
                    
                case 'Mid20RS'
                    keyboard
                    [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    UseCells   = iRS(iSUdps(FirstCell+(0:19)));
                case 'Mid10RS'
                    keyboard
                    [dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    UseCells   = iRS(iSUdps(FirstCell+(0:9)));
                    SUdps      = dps(FirstCell+(0:9));
                    
                case 'BestRS'
                    [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    UseCells   = iRS(iSUdps(1:FirstCell));
                case 'SkipBestRS'
                    [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    UseCells   = iRS(iSUdps((FirstCell+1):end));
                    
                case 'LoudestRS'
                    nSpkRS     = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
                    [~,iSUFR]  = sort(nSpkRS,'descend');
                    UseCells   = iRS(iSUFR(1:FirstCell));
                case 'SkipLoudestRS'
                    nSpkRS     = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
                    [~,iSUFR]  = sort(nSpkRS,'descend');
                    UseCells   = iRS(iSUFR((FirstCell+1):end));
                case 'LowFFRS'
                    keyboard
                    [~,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
                    UseCells   = iRS(iSUdps((FirstCell+1):end));
            end
            
            if exclNonSig
                if numel(UnSig) ~= numel(theseCells), keyboard, end
                UseCells = UseCells(UnSig(iRS));
            end
            if exclSpec
                switch whichStim
                    case 'AC'
                        UseCells = UseCells(CReach(UseCells,:).dprime<1);
                    case 'Speech'
%                         [~,iSUdps] = sort(CReach(UseCells,:).dprime,'descend');
%                         UseCells = UseCells(iSUdps(8:end));
                        UseCells = UseCells(CReach(UseCells,:).dprime<1.4);
                end
            end
            SUdps = CReach(UseCells,:).dprime
            
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if strcmp(whichClass,'OnlyTemp')
                keyboard
            elseif strcmp(whichClass,'Nspk')
                S_AssMat = runSVMclassFR( CTTS(UseCells,:,:,:), BootstrapN, nStim, Dur, length(UseCells), PickTrials, TrainSize, TestSize, KernelType );
            elseif strcmp(whichClass,'Full')
                
                % Train and test classifier
                
%                 [S_AssMat,~,~] = runSVMclass_SUM( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
%                     BootstrapN, nStim, Dur, length(UseCells), PickTrials,TrainSize, TestSize, KernelType, ShuffOpt ); 
                
                [S_AssMat,~,~] = runSVMclass_notNorm( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
                    BootstrapN, nStim, Dur, length(UseCells), PickTrials,TrainSize, TestSize, KernelType, ShuffOpt ); 
            
            elseif strcmp(whichClass,'Sum')
                [S_AssMat,~,~] = runSVMclass_SUM( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
                    BootstrapN, nStim, Dur, length(UseCells), PickTrials,TrainSize, TestSize, KernelType, ShuffOpt );
                
            elseif strcmp(whichClass,'SumAV')
                [S_AssMat,~,~] = runSVMclass_avSUM( CTTS(UseCells,:,:,:), ETTS(UseCells,:,:,:), ...
                    BootstrapN, nStim, Dur, length(UseCells), PickTrials,TrainSize, TestSize, KernelType, ShuffOpt );
                
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
            
            title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s (N=%i)',...
                muPC,dprime,KernelType,whichStim,whichCells,numel(UseCells)))
%             title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s ',...
%                 muPC,dprime,KernelType,whichStim,whichCells))
            
            
            % Save figure
%             savename = sprintf('Res_v%s-%i_%s_%s_exclNS%i_exclSig%i',whichClass,WinEnds(ii)-WinBeg(ii)+1,whichStim,whichCells,exclNonSig,exclSpec);
%             savename = sprintf('Res_v%s-%i_%s_%s_nc%i_fc%i',whichClass,WinEnds(ii)-WinBeg(ii)+1,whichStim,whichCells,NumCells,FirstCell);
            savename = sprintf('Res_%s_%s_nc%i_win%i',whichStim,whichCells,numel(UseCells),WinBeg(ii));
            
            print(hf(ii),fullfile(figsavedir,savename),'-dpdf')
            
            
            %% Save results to master table
            
            mastertablesavename = sprintf('CR_v%s_%s',whichClass,whichCells);
            thistablesavename   = savename;
            
            CR1 = table;
            CR1.figname  = {savename};
            CR1.Stim     = {whichStim};
            CR1.Cells    = {whichCells};
            CR1.exNonSig = exclNonSig;
            CR1.exclSpec = exclSpec;
            CR1.iC       = FirstCell;
            CR1.nC       = numel(UseCells);
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
end %analysis window 

return
keyboard

[~,ConfIntervals] = bootstrap4significance(CR);

idx = [4 3 1 2];
idx = 1:4;

figure;
plot([1:4; 1:4],ConfIntervals(idx,:)','k-','LineWidth',2)
hold on
plot(1:4,CR.dprime(idx),'.k','MarkerSize',30)
xlim([0 5])
ylim([0 4])
print_eps_kp(gcf,fullfile(figsavedir,'ResultPlot'))


% Plot ranked SU vs subpops
pcr_SUvPools_v2

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
