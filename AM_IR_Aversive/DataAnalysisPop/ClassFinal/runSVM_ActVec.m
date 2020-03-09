function [S_AssMat,TrialResults] = runSVM_ActVec(CTTS,BSN,nSt,Dur,nUn,PickTrials,TrainSize,TestSize,KernelType,ShuffOpt)
% AssMat = runSVM_ActVec(Cell_Time_Trial_Stim,BootstrapN,...
%    nStim,Dur,nUns,TrainSize,TestSize,convwin,KernelType)
%
%  Uses activity vectors instead of projection values.
%
%  KP, 2020-02 

S_AssMat     = nan(nSt,nSt,BSN);
TrialResults = [];

h=waitbar(0); 

for iBS = 1:BSN
    
    waitbar(iBS/BSN,h,'Bootstrapping sets of training and testing trials...')
    
    if ShuffOpt
        CTTS = shuffleSpikeTimes(CTTS);
    end
    
    %% If want all cells to have simultaneous trials, define all sets here
    
    if strcmp(PickTrials,'sim')
        %======================
        % SIMULTANEOUS trials
        %======================
        nTrials   = nan(nSt,1);
        Train_trs = nan(nSt,TrainSize);
        Test_trs  = nan(nSt,TestSize);
        for itmp = 1:nSt
            if ~isempty(find(~isnan(mean(CTTS(1,:,:,itmp),2)),1,'last'))
                nTrials(itmp) = find(~isnan(mean(CTTS(1,:,:,itmp),2)),1,'last');
            else
                keyboard
            end
            
            % TRAINING
            Train_trs(itmp,:)  = randperm(nTrials(itmp),TrainSize);
            
            % TESTING
            SetRest = 1:nTrials(itmp);
            SetRest(Train_trs(itmp,:)) = [];
            Test_tr_ids  = randperm(length(SetRest),TestSize);
            Test_trs(itmp,:) = SetRest(Test_tr_ids);
        end
    end
    
    
    %% Preallocate
    
    % Data -- 40 x 2056
    %  observations: trials x stim
    %  predictors: templates x cells
    
    DataTrain   = [];
    DataTest    = [];
    EnvTrain    = [];
    EnvTest     = [];
    TrueTrain   = [];
    TrueTest    = [];
    
    
    %~~~~~~~~~~~~~~~~~~~~
    % For each cell...
    %~~~~~~~~~~~~~~~~~~~~
    for iUn = 1:nUn
        
        Un_DataTrain   = [];
        Un_DataTest    = [];
        Un_EnvTrain    = [];
        Un_EnvTest     = [];
        
        %% Collect PSTHs for stimulus templates
        
        if strcmp(PickTrials,'rand')
            %=================
            %  RANDOM trials
            %=================
            nTrials   = nan(nSt,1);
            Train_trs = nan(nSt,TrainSize-1);
            Test_trs  = nan(nSt,TestSize+1);
            for itmp = 1:nSt
                if ~isempty(find(~isnan(mean(CTTS(iUn,:,:,itmp),2)),1,'last'))
                    nTrials(itmp) = find(~isnan(mean(CTTS(iUn,:,:,itmp),2)),1,'last');
                else
                    keyboard
                end
                
                %=====TRAINING=====
                Train_trs(itmp,:)  = randperm(nTrials(itmp),TrainSize-1);
                
                
                %=====TESTING=====
                SetRest = 1:nTrials(itmp);
                SetRest(Train_trs(itmp,:)) = [];
                Test_tr_ids  = randperm(length(SetRest),TestSize+1);
                Test_trs(itmp,:) = SetRest(Test_tr_ids);
                
            end
        end
        
        
        %% ###############################################################
        %                     Prepare input data
        %  ###############################################################
        
        for ist = 1:nSt   %true stimulus of test trials
            
            %==================
            %  Training data
            %==================
            for itr = Train_trs(ist,:)
                % Get this trial data
                Un_DataTrain  = [ Un_DataTrain; CTTS(iUn,:,itr,ist) ];
                TrueTrain     = [ TrueTrain;    ist ];
            end
            
            %==================
            %    Test data
            %==================
            for itr = Test_trs(ist,:)
                % Get trial data
                Un_DataTest  = [ Un_DataTest; CTTS(iUn,:,itr,ist) ];
                TrueTest     = [ TrueTest;    ist ];
            end
            
        end %ist
        
        if any(any(isnan(Un_DataTrain)))
            keyboard
        end
        DataTrain   = [ DataTrain Un_DataTrain];
        DataTest    = [ DataTest   Un_DataTest];
        
    end %iUn
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %####  NEURAL
    
    % Finish preparing input data
    DataTrain = single(round(DataTrain,2));
    DataTest  = single(round(DataTest,2));
    
    TrueTrain = single(TrueTrain(1:size(DataTrain,1)));
    TrueTest  = single(TrueTest(1:size(DataTest,1)));
    
    InputData = [DataTrain TrueTrain];
    
    % - - - - - - - - - 
    % Train classifier 
    trainedClass = trainSVMClass(InputData,KernelType);
    
    % to get weights (for each learner, one-vs-all)
    % trainedClass.ClassificationSVM.BinaryLearners{1}.Beta
    
    % - - - - - - - - - - - - - 
    % Test on held out data
    Test_FitStim  = trainedClass.predictFcn(DataTest);
    
    confMat = confusionmat(single(TrueTest), Test_FitStim);
    S_AssMat(:,:,iBS) = confMat./sum(confMat,2);
    
    % - - - - - - - - - - - - - 
    % Save stimulus, trial number, and outcome
%     TrialResults = [TrialResults; [TrueTest Test_FitStim Test_trs] ];
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
end %iBS

fprintf('Overall mean PC = %0.1f\n',mean(diag(mean(S_AssMat,3)))*100)
close(h)

end