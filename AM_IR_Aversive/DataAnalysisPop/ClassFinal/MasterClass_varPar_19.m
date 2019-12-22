function MasterClass_varPar_19
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

% close all

varPar       = 'AnDur';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'all'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = [50 75 100 125 150 200 300 400 500];
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrialss  = {'rand'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichIrr     = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BootstrapN   = 600;
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
rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'Cell_Time_Trial_Stim')); 
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
% Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;
if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
    Un_Indices = 1:257;
end

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


%% Prepare to parse data

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end


switch whichIrr
    case 'AC'
        theseStim  = 1:8;
    case 'DB'
        theseStim  = [1:6 9:11];
    case 'Speech'
        theseStim  = 1:size(Cell_Time_Trial_Stim,4);
%         MinStim  = find(min(nTrialMat,[],1)>=minTrs);
end

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>3);



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%##########################################################################
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ii = 1:numel(WinEnds)
    
    AnWin = WinBeg(ii):WinEnds(ii);
    
    for iTrSh = 1:numel(PickTrialss)
        
        PickTrials = PickTrialss{iTrSh};
        
        % Define cells and stimuli
        [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
            whichIrr, whichCells, nTrialMat, UnitData,...
            theseStim, iRS, iNS, minTrs, convwin, AnWin );
        
        
        % Train and test classifier
        AssMat = runSVMclassFR( CTTS, BootstrapN, nStim, Dur, nUns, ...
            PickTrials,TrainSize, TestSize, KernelType );
        
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        %##########################################################################
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        %%
        % Figure settings
        set(0,'DefaultTextInterpreter','none')
        set(0,'DefaultAxesFontSize',18)
        
        % scrsz = get(0,'ScreenSize');     %[left bottom width height]
        % fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
        % shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];
        
        
        %% Plot NEURAL 
        ConfMat = mean(AssMat,3);
        muPC    = mean(diag(ConfMat))*100;
        dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
        
        
        % Plot
        hf(ii) = figure;
        imagesc(ConfMat)
        axis square
        caxis([0 1])
        cmocean('ice') %ice
        % cmocean('curl','pivot',0)
        colorbar
        ylabel('True stim')
        xlabel('Assigned')
        set(gca,'tickdir','out','xtick',1:nStim,'ytick',1:nStim)
        
        title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %s (N=%i)',...
            muPC,dprime,KernelType,whichIrr,whichCells,nUns))
        
        
        % Save figure
        savedir = fullfile(fn.figs,'ClassAM');
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        
%         savename = sprintf('Res_v%s_Train%i_%s_%s_%s',varPar,...
%             TrainSize,whichIrr,whichCells,PickTrials);
        savename = sprintf('Res_v%s-%i_Train%i_conv%i_%s_%s',varPar,WinEnds(ii)-WinBeg(ii)+1,TrainSize,tau,whichIrr,whichCells);
        
        print(hf(ii),fullfile(savedir,savename),'-dpdf')
        
        
        %% Save results to master table
        
        mastertablesavename = sprintf('CR_v%s_%s_%s',varPar,whichIrr,whichCells);
%         mastertablesavename = sprintf('CR_v%s_%s',varPar,whichIrr);
        thistablesavename   = savename;
        
        CR1 = table;
        CR1.figname  = {savename};
        CR1.Stim     = {whichIrr};
        CR1.Cells    = {whichCells};
        CR1.iC       = 0;
        CR1.Results  = {AssMat};
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
            save(fullfile(savedir,thistablesavename),'CR1','-v7.3')
            
        else % Save updated table
            
            % Concatenate new data
            CR = q.CR;
            CR = [CR; CR1];
            
            save(fullfile(savedir,mastertablesavename),'CR','-v7.3')
            save(fullfile(savedir,thistablesavename),'CR1','-v7.3')
        end
        
        
    end % iTrSh
end %vary classification parameter


keyboard

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


% 
figure; 
plot([CR.WinEnd]-500,[CR.dprime],'ok')
xlabel('Time from onset (ms)')
ylabel('d''')
ylim([0 3])
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
