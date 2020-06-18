function [S_AssMat,E_AssMat,TrialResults] = runSVMclass_notNorm(CTTS,ETTS,BSN,nSt,Dur,nUn,PickTrials,TrainSize,TestSize,KernelType,ShuffOpt)
% AssMat = runSVMclass(Cell_Time_Trial_Stim,BootstrapN,...
%    nStim,Dur,nUns,TrainSize,TestSize,convwin,KernelType)
%
%  Called by MasterClass
% 

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
            Train_trs = nan(nSt,TrainSize);
            Test_trs  = nan(nSt,TestSize);
            for itmp = 1:nSt
                if ~isempty(find(~isnan(mean(CTTS(iUn,:,:,itmp),2)),1,'last'))
                    nTrials(itmp) = find(~isnan(mean(CTTS(iUn,:,:,itmp),2)),1,'last');
                else
                    keyboard
                end
                
                %=====TRAINING=====
                Train_trs(itmp,:)  = randperm(nTrials(itmp),TrainSize);
                
                
                %=====TESTING=====
                SetRest = 1:nTrials(itmp);
                SetRest(Train_trs(itmp,:)) = [];
                Test_tr_ids  = randperm(length(SetRest),TestSize);
                Test_trs(itmp,:) = SetRest(Test_tr_ids);
                
            end
        end
        
% % %         CTTS = shuffleSpikeTimes(CTTS);
        
        %=================
        %   TEMPLATES
        %=================
        T  = nan(nSt,Dur);
        TE = nan(nSt,Dur);
        for itmp = 1:nSt
            % Get template for comparison to test trial (with 1 random
            % trial left out, to match ntrials in the PSTH template for
            % computing training data corr values) 
            T(itmp,:)   = mean(CTTS(iUn,:,Train_trs(itmp,randperm(TrainSize,TrainSize-1)),itmp),3);
            TE(itmp,:)  = mean(ETTS(iUn,:,Train_trs(itmp,randperm(TrainSize,TrainSize-1)),itmp),3);
        end 
        
        
        %% ###############################################################
        %                     Prepare input data
        %  ###############################################################
        
        for ist = 1:nSt   %true stimulus of test trials
            
            %==================
            %  Training data
            %==================
            for itr = Train_trs(ist,:)
                
                % Get trial data
                trial_data = CTTS(iUn,:,itr,ist);
                trEnv_data = ETTS(iUn,:,itr,ist);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . 
                % . Dot product of this test trial with each template
                R  = nan(1,nSt);
                RE = nan(1,nSt);
                for itmp = 1:nSt
                    
                    % Get template data
%                               Nm1_PSTH = mean( CTTS(iUn,:,Train_trs(itmp,Train_trs(itmp,:)~=itr),itmp) ,3);
%                               Nm1_PETH = mean( ETTS(iUn,:,Train_trs(itmp,Train_trs(itmp,:)~=itr),itmp) ,3);
                    Nm1_PSTH = mean( CTTS(iUn,:,Train_trs(itmp,Train_trs(ist,:)~=itr),itmp) ,3);
                    Nm1_PETH = mean( ETTS(iUn,:,Train_trs(itmp,Train_trs(ist,:)~=itr),itmp) ,3);
                    
                    % Calculate dot product
                    if sum(trial_data)>0 && sum(Nm1_PSTH)>0 
                        R(itmp) = (trial_data * Nm1_PSTH') ;
                    else
                        R(itmp) = 0;
                    end
                    RE(itmp) = (trEnv_data * Nm1_PETH') ;
                end
                Un_DataTrain  = [ Un_DataTrain;   R ];
                Un_EnvTrain   = [ Un_EnvTrain;   RE ];
                TrueTrain     = [ TrueTrain;    ist ];
            end
            
            %==================
            %    Test data
            %==================
            for itr = Test_trs(ist,:)
                
                % Get trial data
                trial_data = CTTS(iUn,:,itr,ist);
                trEnv_data = ETTS(iUn,:,itr,ist);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . 
                % . Dot product of this test trial with each template
                R  = nan(1,nSt);
                RE = nan(1,nSt);
                for itmp = 1:nSt
                    if sum(trial_data)>0 && sum(T(itmp,:))>0 %&& mod(itr,2)>0
                        R(itmp) = (trial_data * T(itmp,:)') ;
                    else
                        R(itmp) = 0;
                    end
                    RE(itmp) = (trEnv_data * TE(itmp,:)') ;
                end
                
                Un_DataTest  = [ Un_DataTest;   R ];
                Un_EnvTest   = [ Un_EnvTest;   RE ];
                TrueTest     = [ TrueTest;    ist ];
            end
            
        end %ist
        
        if any(any(isnan(Un_DataTrain)))
            keyboard
        end
        DataTrain   = [ DataTrain Un_DataTrain];
        DataTest    = [ DataTest   Un_DataTest];
        EnvTrain    = [ EnvTrain   Un_EnvTrain];
        EnvTest     = [ EnvTest     Un_EnvTest];
        
    end %iUn
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %####  NEURAL
    
    % Finish preparing input data
    DataTrain = single(round(DataTrain,4));
    DataTest  = single(round(DataTest,4));
    
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
    TrialResults = [TrialResults; [TrueTest Test_FitStim Test_trs] ];
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %####  ENVELOPE  
    
    % Finish preparing input data
% % %     InputData = [];
% % %     InputData = [EnvTrain TrueTrain];
% % %     
% % %     % - - - - - - - - - 
% % %     % Train classifier 
% % %     trainedClass = trainSVMClass(InputData,KernelType);
% % %     
% % %     % - - - - - - - - - - - - - 
% % %     % Test on held out data
% % %     Test_FitStim  = trainedClass.predictFcn(EnvTest);
% % %     
% % %     confMat = confusionmat(TrueTest, Test_FitStim);
% % %     E_AssMat(:,:,iBS) = confMat./sum(confMat,2);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
end %iBS

fprintf('Overall mean PC = %0.1f, Env %0.1f%%\n',mean(diag(mean(S_AssMat,3)))*100,mean(diag(mean(E_AssMat,3)))*100)
close(h)

end