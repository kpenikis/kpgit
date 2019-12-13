function AssMat = runSVMclass(CTTS,BSN,nSt,Dur,nUn,PickTrials,TrainSize,TestSize,KernelType)
% AssMat = runSVMclass(Cell_Time_Trial_Stim,BootstrapN,...
%    nStim,Dur,nUns,TrainSize,TestSize,convwin,KernelType)
%
%  Called by MasterClass
% 

AssMat = nan(nSt,nSt,BSN);

h=waitbar(0);

for iBS = 1:BSN
    
    waitbar(iBS/BSN,h,'Bootstrapping sets of training and testing trials...')
    
    
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
    TrueTrain   = [];
    TrueTest    = [];
    
    
    %~~~~~~~~~~~~~~~~~~~~
    % For each cell...
    %~~~~~~~~~~~~~~~~~~~~
    for iUn = 1:nUn
        
        Un_DataTrain   = [];
        Un_DataTest    = [];
        
        
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
                
                % TRAINING
                Train_trs(itmp,:)  = randperm(nTrials(itmp),TrainSize);
                
                % TESTING
                SetRest = 1:nTrials(itmp);
                SetRest(Train_trs(itmp,:)) = [];
                Test_tr_ids  = randperm(length(SetRest),TestSize);
                Test_trs(itmp,:) = SetRest(Test_tr_ids);
            end
        end
        
        
        %=================
        %   TEMPLATES
        %=================
        T = nan(nSt,Dur);
        for itmp = 1:nSt
            
            % N-1 trials
            PSTH_trs = 1:nTrials(itmp);
            PSTH_trs([Train_trs(itmp,:) Test_trs(itmp,:)]) = [];
            
            % Make template
            T(itmp,:)  = mean(CTTS(iUn,:,PSTH_trs,itmp),3);
            
        end %itmp
        
        
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
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . 
                % . Dot product of this test trial with each template
                R = nan(1,size(T,1));
                for itmp = 1:nSt
                    if sum(trial_data)>0 && sum(T(itmp,:))>0 %&& mod(itr,2)>0
                        R(itmp) = (trial_data*T(itmp,:)')/(norm(trial_data)*norm(T(itmp,:)));
                    else
                        R(itmp) = 0;
                    end
                end
                Un_DataTrain  = [ Un_DataTrain;   R ];
                TrueTrain  = [ TrueTrain; ist ];
            end
            
            %==================
            %    Test data
            %==================
            for itr = Test_trs(ist,:)
                
                % Get trial data
                trial_data = CTTS(iUn,:,itr,ist);
                
                % . .  . . . . . . . . . . . . . . . . . . . . . . . . . . 
                % . Dot product of this test trial with each template
                R = nan(1,size(T,1));
                for itmp = 1:nSt
                    if sum(trial_data)>0 && sum(T(itmp,:))>0 %&& mod(itr,2)>0
                        R(itmp) = (trial_data*T(itmp,:)')/(norm(trial_data)*norm(T(itmp,:)));
                    else
                        R(itmp) = 0;
                    end
                end
                
                Un_DataTest  = [ Un_DataTest;   R ];
                TrueTest  = [ TrueTest; ist ];
            end
            
        end %ist
        
        if any(any(isnan(Un_DataTrain)))
            keyboard
        end
        DataTrain   = [ DataTrain Un_DataTrain];
        DataTest    = [ DataTest  Un_DataTest];
        
    end %iUn
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Finish preparing input data
    TrueTrain = TrueTrain(1:size(DataTrain,1));
    TrueTest  = TrueTest(1:size(DataTest,1));
    
    InputData = [DataTrain TrueTrain];
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Train classifier
    trainedClass = trainSVMClass(InputData,KernelType);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Test on held out data
    Test_FitStim  = trainedClass.predictFcn(DataTest);
    
    confMat = confusionmat(TrueTest, Test_FitStim);
    AssMat(:,:,iBS) = confMat./sum(confMat,2);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
end %iBS

fprintf('Overall mean PC = %0.1f%%, %i tr Training\n',mean(diag(mean(AssMat,3)))*100,TrainSize)
close(h)

end