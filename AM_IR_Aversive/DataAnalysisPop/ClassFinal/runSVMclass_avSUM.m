function [S_AssMat,E_AssMat,TrialResults] = runSVMclass_avSUM(CTTS,ETTS,BSN,nSt,Dur,nUn,PickTrials,TrainSize,TestSize,KernelType,ShuffOpt)
% AssMat = runSVMclass_avSUM
%
%  Called by MasterClass programs
%
%  This version uses summed activity vector as input data. 
%  runSVMclass_SUM uses projection values.
%
%  KP, 2020-07 

S_AssMat     = nan(nSt,nSt,BSN);
E_AssMat     = [];
TrialResults = [];

h=waitbar(0);

for iBS = 1:BSN
    
    waitbar(iBS/BSN,h,'Bootstrapping sets of training and testing trials...')
    
    if ShuffOpt
        CTTS = shuffleSpikeTimes(CTTS);
    end
    
    %% If want all cells to have simultaneous trials, define all sets here
                    
     FiltDataTrain = nan(nUn,size(CTTS,2),TrainSize-1,nSt);
     FiltDataTest  = nan(nUn,size(CTTS,2),TestSize+1,nSt);
        
    if strcmp(PickTrials,'sim')
        keyboard
        %======================
        % SIMULTANEOUS trials
        %======================
        nTrials   = nan(nSt,1);
        Train_trs = nan(nSt,TrainSize-1);
        Test_trs  = nan(nSt,TestSize+1);
        for ist = 1:nSt
            if ~isempty(find(~isnan(mean(CTTS(1,:,:,ist),2)),1,'last'))
                nTrials(ist) = find(~isnan(mean(CTTS(1,:,:,ist),2)),1,'last');
            else
                keyboard
            end
            
            %=====TRAINING=====
            Train_trs(ist,:)  = randperm(nTrials(ist),TrainSize-1);
            
            FiltDataTrain(:,:,:,ist) = CTTS(:,:,Train_trs(ist,:),ist);
            Train_trs(ist,:)  = 1:size(Train_trs,2);
            
            %=====TESTING=====
            SetRest = 1:nTrials(ist);
            SetRest(Train_trs(ist,:)) = [];
            Test_tr_ids  = randperm(length(SetRest),TestSize+1);
            Test_trs(ist,:) = SetRest(Test_tr_ids);
            
            FiltDataTest(iUn,:,:,ist) = CTTS(:,:,Test_trs(ist,:),ist);
            Test_trs(ist,:)  = 1:size(Test_trs,2);
            
        end
        
    elseif strcmp(PickTrials,'rand')
        
        %=================
        %  RANDOM trials
        %=================
        for iUn = 1:nUn
            
            nTrials   = nan(nSt,1);
            Train_trs = nan(nSt,TrainSize-1);
            Test_trs  = nan(nSt,TestSize+1);
            for ist = 1:nSt
                if ~isempty(find(~isnan(mean(CTTS(iUn,:,:,ist),2)),1,'last'))
                    nTrials(ist) = find(~isnan(mean(CTTS(iUn,:,:,ist),2)),1,'last');
                else
                    keyboard
                end
                
                %=====TRAINING=====
                Train_trs(ist,:)  = randperm(nTrials(ist),TrainSize-1);
                
                if any(any(isnan(CTTS(iUn,:,Train_trs(ist,:),ist))))
                    keyboard
                end
                FiltDataTrain(iUn,:,:,ist) = CTTS(iUn,:,Train_trs(ist,:),ist);
%                 Train_trs(ist,:)  = 1:size(Train_trs,2);
                
                
                %=====TESTING=====
                SetRest = 1:nTrials(ist);
                SetRest(Train_trs(ist,:)) = [];
                Test_tr_ids  = randperm(length(SetRest),TestSize+1);
                Test_trs(ist,:) = SetRest(Test_tr_ids);
                
                FiltDataTest(iUn,:,:,ist) = CTTS(iUn,:,Test_trs(ist,:),ist);
%                 Test_trs(ist,:)  = 1:size(Test_trs,2);
                
            end %ist
            
        end %iUn (rand)
        
    end %sim or rand
    
    for ist = 1:nSt
        Train_trs(ist,:) = 1:size(Train_trs,2);
        Test_trs(ist,:)  = 1:size(Test_trs,2);
    end
    
    
    % Sum over cells !!
    
    FiltDataTrain = sum(FiltDataTrain,1);
    FiltDataTest  = sum(FiltDataTest,1);
    
    
    %% ###############################################################
    %                     Prepare input data
    %  ###############################################################
    
    % Preallocate
%     Un_DataTrain   = nan(size(FiltDataTrain,3)*size(FiltDataTrain,4),size(FiltDataTrain,2)); % Un_DataTrain [ Ntr*Nst time ]
%     Un_DataTest    = nan(size(FiltDataTest,3)*size(FiltDataTest,4),size(FiltDataTest,2));
    
    Un_DataTrain = [];
    Un_DataTest  = [];
    DataTrain    = [];
    DataTest     = [];
    TrueTrain    = [];
    TrueTest     = [];
    
    for ist = 1:nSt   %true stimulus of test trials
        
        %==================
        %  Training data
        %==================
        for itr = Train_trs(ist,:)
            % Get this trial data
            Un_DataTrain  = [ Un_DataTrain; FiltDataTrain(:,:,itr,ist) ];
            TrueTrain     = [ TrueTrain;    ist ];
        end
        
        %==================
        %    Test data
        %==================
        for itr = Test_trs(ist,:)
            % Get trial data
            Un_DataTest  = [ Un_DataTest; FiltDataTest(:,:,itr,ist) ];
            TrueTest     = [ TrueTest;    ist ];
        end
        
    end %ist
    
    if any(any(isnan(Un_DataTrain)))
        keyboard
    end
    DataTrain   = [ DataTrain Un_DataTrain];
    DataTest    = [ DataTest   Un_DataTest];
    
    
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

fprintf('Overall mean PC = %0.1f%%\n',mean(diag(mean(S_AssMat,3)))*100)
close(h)

end