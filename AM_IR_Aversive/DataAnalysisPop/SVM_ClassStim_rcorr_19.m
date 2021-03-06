function SVM_ClassStim_rcorr_19
% SVM_ClassStim_rcorr
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


% must be the wrong organization of input data

% 1st -- 6920 x 8
%  observations: trials x stim x cells 
%  predictors: templates

% Now do -- 40 x 2056
%  observations: trials x stim
%  predictors: templates x cells 


% close all

whichIrr    = 'AC';

whichCells  = 'each'; 
PickTrials  = 'rand';

Un_Index    = [40 115];

BootstrapN  = 500;
KernelType  = 'linear';



tau     = 10;
lambda  = 1/tau;
winlen  = 500;
convwin = exp(-lambda*(1:winlen));

PSTHsize  = 8;
TrainSize = 5;
TestSize  = 5;
minTrs    = PSTHsize + TrainSize + TestSize;


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(savedir,'Cell_Time_Trial_Stim_simtrs'));
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
    Un_Indices = 1:257;
end

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>3);


%% Prepare input data

% Can't do entire population (or complete groups) with ALL stimuli
% including Irr. Can do AC and DB separately.

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end


%% Define cells and stimuli

ACstim  = 1:8;
DBstim  = [1:6 9:11];

switch whichCells
    case {'all' 'RSNS' 'each'}                                         % All Cells
        switch whichIrr
            case 'AC'
                theseCells = find(all(nTrialMat(:,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = find(all(nTrialMat(:,DBstim)>=minTrs,2));
        end
        
        if ~strcmp(whichCells,'all')
            iiRS = intersect(theseCells,iRS);
            iiNS = intersect(theseCells,iNS);
            theseCells = sort([iiRS; iiNS]);
            
            [~,RS_CTTS] = intersect(theseCells,iiRS);
            [~,NS_CTTS] = intersect(theseCells,iiNS);
        end
        
        
    case {'Mar28-AM' 'Mar30-AM' 'Apr02-AM' 'Apr11-AM'}            % Session
        
        switch whichIrr
            case 'AC'
                theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,ACstim)>=minTrs,2) );
            case 'DB'
                theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,DBstim)>=minTrs,2) );
        end
        
        
    case 'select'                                  % One or subset of cells
        keyboard
        theseCells = Un_Index;
        
        
    case 'RS'                                                          % RS
        switch whichIrr
            case 'AC'
                theseCells = iRS(all(nTrialMat(iRS,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = iRS(all(nTrialMat(iRS,DBstim)>=minTrs,2));
        end
        
        
    case 'NS'                                                          % NS
        switch whichIrr
            case 'AC'
                theseCells = iNS(all(nTrialMat(iNS,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = iNS(all(nTrialMat(iNS,DBstim)>=minTrs,2));
        end
end


switch whichIrr
    case 'AC'
        Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(theseCells,:,:,ACstim);
        nTrialMat = nTrialMat(theseCells,ACstim);
    case 'DB'
        Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(theseCells,:,:,DBstim);
        nTrialMat = nTrialMat(theseCells,DBstim);
end

% Complete
nUns   = size(Cell_Time_Trial_Stim,1);
Dur    = size(Cell_Time_Trial_Stim,2);
minTrStim = min(nTrialMat,[],1);
nStim  = size(Cell_Time_Trial_Stim,4);


%% Set size of training and testing data sets
% 
% TrainSize = 5;
% TestSize  = 5;
% TestSize  = min(8,min(minTrStim)-TrainSize);


%%

CRall = table;


for iUn = 1:nUns

fprintf('Classifying sess %s, clu %i\n',UnitData(theseCells(iUn)).Session,UnitData(theseCells(iUn)).Clu)

AssignmentMat = nan(nStim,nStim,BootstrapN);

h=waitbar(0);

for iBS = 1:BootstrapN
    
    waitbar(iBS/BootstrapN,h,'Bootstrapping sets of training and testing trials...')
    
    
    %% If want all cells to have simultaneous trials, define all sets here
    
    if strcmp(PickTrials,'sim')
        
        nTrials   = nan(nStim,1);
        PSTH_trs  = nan(nStim,PSTHsize);
        
        Train_trs = nan(nStim,TrainSize);
        Test_trs  = nan(nStim,TestSize);
        
        % TEMPLATES
        for itmp = 1:nStim
            if ~isempty(find(~isnan(mean(Cell_Time_Trial_Stim(1,:,:,itmp),2)),1,'last'))
                nTrials(itmp) = find(~isnan(mean(Cell_Time_Trial_Stim(1,:,:,itmp),2)),1,'last');
            else
                keyboard
            end
            PSTH_trs(itmp,:)  = randperm(nTrials(itmp),PSTHsize);
            
            % TRAINING
            SetRest = 1:nTrials(itmp);
            SetRest(PSTH_trs(itmp,:)) = [];
            Train_tr_ids  = randperm(length(SetRest),TrainSize);
            Train_trs(itmp,:) = SetRest(Train_tr_ids);
            if any(intersect(Train_trs(itmp,:),PSTH_trs(itmp,:)))
                keyboard
            end
            
            % TESTING
            SetRest = 1:nTrials(itmp);
            SetRest([PSTH_trs(itmp,:) Train_trs(itmp,:)]) = [];
            Test_tr_ids  = randperm(length(SetRest),TestSize);
            Test_trs(itmp,:) = SetRest(Test_tr_ids);
            if any(intersect(Test_trs(itmp,:),[PSTH_trs(itmp,:) Train_trs(itmp,:)]))
                keyboard
            end
            
        end
    end
    
    
    %% Preallocate
    
    % Data -- 40 x 2056
    %  observations: trials x stim
    %  predictors: templates x cells
    
    TrainingData   = [];
    TestingData    = [];
    Train_TrueStim = [];
    Test_TrueStim  = [];
    
    
    %~~~~~~~~~~~~~~~~~~~~~~
    % For each cell...
    %~~~~~~~~~~~~~~~~~~~~~~
%     for iUn = 1:nUns
        
        UnTrainingData   = [];
        UnTestingData    = [];
        
        
        %% Collect PSTHs for stimulus templates
        
        if strcmp(PickTrials,'rand')
            
            nTrials   = nan(nStim,1);
            PSTH_trs  = nan(nStim,PSTHsize);
            T         = nan(nStim,Dur);
            
            Train_trs = nan(nStim,TrainSize);
            Test_trs  = nan(nStim,TestSize);
            
            % TEMPLATES
            for itmp = 1:nStim
                if ~isempty(find(~isnan(mean(Cell_Time_Trial_Stim(iUn,:,:,itmp),2)),1,'last'))
                    nTrials(itmp) = find(~isnan(mean(Cell_Time_Trial_Stim(iUn,:,:,itmp),2)),1,'last');
                else
                    keyboard
                end
                PSTH_trs(itmp,:)  = randperm(nTrials(itmp),PSTHsize);
                
                %-----convolved templates
                conv_data = [];
                conv_data = conv(mean(Cell_Time_Trial_Stim(iUn,:,PSTH_trs(itmp,:),itmp),3),convwin);
                T(itmp,:)  = conv_data(1:Dur);
                
                
                % TRAINING
                SetRest = 1:nTrials(itmp);
                SetRest(PSTH_trs(itmp,:)) = [];
                Train_tr_ids  = randperm(length(SetRest),TrainSize);
                Train_trs(itmp,:) = SetRest(Train_tr_ids);
                if any(intersect(Train_trs(itmp,:),PSTH_trs(itmp,:)))
                    keyboard
                end
                
                
                % TESTING
                SetRest = 1:nTrials(itmp);
                SetRest([PSTH_trs(itmp,:) Train_trs(itmp,:)]) = [];
                Test_tr_ids  = randperm(length(SetRest),TestSize);
                Test_trs(itmp,:) = SetRest(Test_tr_ids);
                if any(intersect(Test_trs(itmp,:),[PSTH_trs(itmp,:) Train_trs(itmp,:)]))
                    keyboard
                end
                
            end %itmp
            
            
        else % just creat Templates for this Cell
            
            T = nan(nStim,Dur);
            
            % TEMPLATES
            for itmp = 1:nStim
                %-----convolved templates
                conv_data = [];
                conv_data = conv(mean(Cell_Time_Trial_Stim(iUn,:,PSTH_trs(itmp,:),itmp),3),convwin);
                T(itmp,:)  = conv_data(1:Dur);
            end %itmp
            
        end
        
%         nTrials   = nan(nStim,1);
%         PSTH_trs  = nan(nStim,PSTHsize);
%         T         = nan(nStim,Dur);
%         
%         for itmp = 1:nStim
%             if ~isempty(find(~isnan(mean(Cell_Time_Trial_Stim(iUn,:,:,itmp),2)),1,'last'))
%                 nTrials(itmp) = find(~isnan(mean(Cell_Time_Trial_Stim(iUn,:,:,itmp),2)),1,'last');
%             else
%                 keyboard
%             end
%             
%             PSTH_trs(itmp,:)  = randperm(nTrials(itmp),PSTHsize);
%             
%             %-----convolved templates
%             conv_data = [];
%             conv_data = conv(mean(Cell_Time_Trial_Stim(iUn,:,PSTH_trs(itmp,:),itmp),3),convwin);
%             T(itmp,:)  = conv_data(1:Dur);
%         end
        
        
        
        %% Now gather and calculate input data
        
        for ist = 1:nStim   %true stimulus of test trials
            
            % Collect trials for training set
%             SetRest = 1:nTrials(ist);
%             SetRest(PSTH_trs(ist,:)) = [];
%             Train_tr_ids  = randperm(length(SetRest),TrainSize);
%             Train_trs = SetRest(Train_tr_ids);
%             if any(intersect(Train_trs,PSTH_trs(ist,:)))
%                 keyboard
%             end
            
            for itr = Train_trs(ist,:)
                
                conv_data = [];
                
                %-----1ms resolution
                %             add_data = mean(Cell_Time_Trial_Stim(:,:,itr,ist),1);
                
                %-----convolved
                conv_data = conv(Cell_Time_Trial_Stim(iUn,:,itr,ist),convwin);
                conv_data = conv_data(1:Dur);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % . Dot product of this test trial with each template
                R = nan(1,size(T,1));
                for itmp = 1:nStim
                    if sum(conv_data)>0 && sum(T(itmp,:))>0
                        R(itmp) = (conv_data*T(itmp,:)')/(norm(conv_data)*norm(T(itmp,:)));
                    else
                        R(itmp) = 0;
                    end
                end
                UnTrainingData  = [ UnTrainingData;   R ];
                Train_TrueStim  = [ Train_TrueStim; ist ];
                
            end
            
            % Collect trials for test set
%             SetRest = 1:nTrials(ist);
%             SetRest([PSTH_trs(ist,:) Train_trs]) = [];
%             Test_tr_ids  = randperm(length(SetRest),TestSize);
%             Test_trs = SetRest(Test_tr_ids);
%             if any(intersect(Test_trs,[PSTH_trs(ist,:) Train_trs]))
%                 keyboard
%             end
            
            
            for itr = Test_trs(ist,:)
                conv_data = [];
                
                %-----1ms resolution
                %             add_data = mean(Cell_Time_Trial_Stim(:,:,itr,ist),1);
                
                %-----convolved
                conv_data = conv(Cell_Time_Trial_Stim(iUn,:,itr,ist),convwin);
                conv_data = conv_data(1:Dur);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % . Dot product of this test trial with each template
                R = nan(1,size(T,1));
                for itmp = 1:nStim
                    if sum(conv_data)>0 && sum(T(itmp,:))>0
                        R(itmp) = (conv_data*T(itmp,:)')/(norm(conv_data)*norm(T(itmp,:)));
                    else
                        R(itmp) = 0;
                    end
                end
                
                UnTestingData  = [ UnTestingData;   R ];
                Test_TrueStim  = [ Test_TrueStim; ist ];
            end
            
        end %ist
        
        if any(any(isnan(UnTrainingData)))
            keyboard
        end
        TrainingData   = [ TrainingData UnTrainingData];
        TestingData    = [ TestingData  UnTestingData];
        
%     end %iUn
        
    Train_TrueStim = Train_TrueStim(1:40);
    Test_TrueStim  = Test_TrueStim(1:40);
    
    inputData     = [];
    inputData = [TrainingData Train_TrueStim];
    
    
    % - - - - - - - - -
    % Train classifier
    trainedClass = trainSVMClass(inputData,KernelType);
    
    % - - - - - - - - - - - -
    % Test on held out data
    Test_FitStim  = trainedClass.predictFcn(TestingData);
    
    confMat = confusionmat(Test_TrueStim, Test_FitStim);
    AssignmentMat(:,:,iBS) = confMat./sum(confMat,2);
    
end %iBS

fprintf('Overall mean PC = %0.1f%%, %i tr Training\n',mean(diag(mean(AssignmentMat,3)))*100,TrainSize)
close(h)


%%
% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

% scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
% shortscreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/2];


ConfMat = mean(AssignmentMat,3);
muPC = mean(diag(ConfMat))*100;
dprime = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1); 


% Plot
hf = figure;
imagesc(ConfMat)
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
title(sprintf('%0.1f%%, d''=%0.2f\n%s SVM (%s)  |  %i/%i/%i trs  |  %s (N=%i)',muPC,dprime,KernelType,whichIrr,PSTHsize,TrainSize,TestSize,whichCells,nUns))
% title(sprintf('%s SVM  (iUn %s, cat)  PC=%0.1f%%  %i/%i trs, conv %i',KernelType,num2str(Un_Index),mean(diag(mean(AssignmentMat,3)))*100,TrainSize,TestSize,length(convwin)))

% Save figure
savedir = fullfile(fn.figs,'StimClassRcorr');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('Results_Train%i_conv%i_%s_%s_%i',TrainSize,tau,whichIrr,whichCells,theseCells(iUn));
% savename = sprintf('Results_%s_Train%i_conv%i_iUns%icat%i',KernelType,TrainSize,length(convwin),Un_Index);

print(fullfile(savedir,savename),'-dpdf')

close(hf);


%% Save results to master table 

CR1 = table;
CR1.figname  = {savename};
CR1.Stim     = {whichIrr};
CR1.Cells    = {whichCells};
CR1.iC       = theseCells(iUn);
CR1.Results  = {AssignmentMat};
CR1.PC       = muPC;
CR1.dprime   = dprime;
CR1.trials   = {PickTrials};
CR1.rcorrOpt = {'norm'};
CR1.nTemp    = PSTHsize;
CR1.nTrain   = TrainSize;
CR1.nTest    = TestSize;
CR1.conv     = {'exp'};
CR1.tau      = tau;

CRall = [CRall; CR1];

end %iUn


% Load saved table
tablesavename = sprintf('CR_%s_%s',whichIrr,whichCells);
clear q;
try
q=load(fullfile(savedir,tablesavename));
end


if ~exist('q','var')
    
    % Save new table
    CR = CR1;
    save(fullfile(savedir,tablesavename),'CR','-v7.3')
    
else % Save updated table
    
    % Concatenate new data
    CR = q.CR;
    CR = [CR; CR1];
    
    save(fullfile(savedir,tablesavename),'CR','-v7.3')
end

keyboard


end
