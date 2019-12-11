function SVM_ClassStim_tsPop
% SVM_ClassStim_tsPop
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%
%
%  KP, 2019-12
%

% Try 500 iterations.
% Next, try setting standardize to true. 


% close all

whichIrr    = 'AC';
whichCells  = 'RSNScat'; %all select
Un_Index    = [40 115];

BootstrapN = 500;
KernelType = 'linear';

bs_gaus     = 10;
% convwin     = ones(1,convwin).*(1/convwin);
convwin = gausswin(bs_gaus);
convwin = convwin-min(convwin);
convwin = convwin/sum(convwin);

minTrs      = 16;


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
    case {'all' 'RSNScat' 'RSNSsum'}   % All Cells
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
        
        
    case 'select' % One or subset of cells
        keyboard
        theseCells = Un_Index;
        
    case 'RS'
        switch whichIrr
            case 'AC'
                theseCells = iRS(all(nTrialMat(iRS,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = iRS(all(nTrialMat(iRS,DBstim)>=minTrs,2));
        end
        
    case 'NS'
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
minTrStim = min(nTrialMat,[],1);
nStim  = size(Cell_Time_Trial_Stim,4);


%% Set size of training and testing data sets

TrainSize = 15;
TestSize  = min(8,min(minTrStim)-TrainSize);


%%

AssignmentMat = nan(nStim,nStim,BootstrapN);

h=waitbar(0);

for iBS = 1:BootstrapN
    
    waitbar(iBS/BootstrapN,h,'Bootstrapping sets of training and testing trials...')
    
    TrainingData   = [];
    TestingData    = [];
    Train_TrueStim = [];
    Test_TrueStim  = [];
    
    inputData      = [];
    
    for ist = 1:nStim
        
        nTr = find(~isnan(mean(mean(Cell_Time_Trial_Stim(:,:,:,ist),2),1)),1,'last');
        
        if isempty(nTr)
            continue
        end
        
        % Collect training set
        Train_trs = randperm(nTr,TrainSize);
        
        for itr = Train_trs
            
            add_data = [];
            
            if strcmp(whichCells,'RSNScat')
                %concatenate mean of RS and NS cells
                add_data = [mean(Cell_Time_Trial_Stim(NS_CTTS,:,itr,ist),1) mean(Cell_Time_Trial_Stim(RS_CTTS,:,itr,ist),1)];
                add_data = conv(add_data,convwin,'same');
                add_data = add_data(2:2:end);
                
            else
            %-----1ms resolution
%             add_data = mean(Cell_Time_Trial_Stim(:,:,itr,ist),1);
            %-----convolved
            add_data = conv(mean(Cell_Time_Trial_Stim(:,:,itr,ist),1),convwin,'same');
            
            %concatenate cells
%             for iUn = 1:nUns
%                 add_data = [add_data Cell_Time_Trial_Stim(iUn,:,itr,ist)];
%             end
%             add_data = conv(add_data,convwin,'same');
            end
            
            TrainingData   = [ TrainingData; add_data ];
            Train_TrueStim = [ Train_TrueStim; ist ];
        end
        
        % Collect testing set
        TestSet = 1:nTr;
        TestSet(Train_trs) = [];
        Test_tr_ids  = randperm(length(TestSet),TestSize);
        Test_trs = TestSet(Test_tr_ids);
        if any(intersect(Test_trs,Train_trs))
            keyboard
        end
        
        for itr = Test_trs
            add_data = [];
            
            if strcmp(whichCells,'RSNScat')
                %concatenate mean of RS and NS cells
                add_data = [mean(Cell_Time_Trial_Stim(NS_CTTS,:,itr,ist),1) mean(Cell_Time_Trial_Stim(RS_CTTS,:,itr,ist),1)];
                add_data = conv(add_data,convwin,'same');
                add_data = add_data(2:2:end);
                
            else
            %-----1ms resolution
%             add_data = mean(Cell_Time_Trial_Stim(:,:,itr,ist),1);
            %-----convolved
            add_data = conv(mean(Cell_Time_Trial_Stim(:,:,itr,ist),1),convwin,'same');
            
            %concatenate cells
%             for iUn = 1:nUns
%                 add_data = [add_data Cell_Time_Trial_Stim(iUn,:,itr,ist)];
%             end
%             add_data = conv(add_data,convwin,'same');
            end
            
            TestingData   = [ TestingData; add_data ];
            Test_TrueStim = [ Test_TrueStim; ist ];
        end
        
    end %ist
    
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


% Plot
hf = figure;
imagesc(mean(AssignmentMat,3))
axis square
caxis([0 1])
cmocean('ice') %ice
% cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
title(sprintf('%s SVM  (%s, %s)  PC=%0.1f%%  %i/%i trs, N=%i',KernelType,whichIrr,whichCells,mean(diag(mean(AssignmentMat,3)))*100,TrainSize,TestSize,nUns))
% title(sprintf('%s SVM  (iUn %s, cat)  PC=%0.1f%%  %i/%i trs, conv %i',KernelType,num2str(Un_Index),mean(diag(mean(AssignmentMat,3)))*100,TrainSize,TestSize,length(convwin)))

% Save figure
savedir = fullfile(fn.figs,'SVMclass');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('Results_%s_Train%i_conv%i_%s_%s_2ms',KernelType,TrainSize,length(convwin),whichIrr,whichCells);
% savename = sprintf('Results_%s_Train%i_conv%i_iUns%icat%i',KernelType,TrainSize,length(convwin),Un_Index);

print(fullfile(savedir,savename),'-dpdf')


keyboard




%%

% save data too ?
save(fullfile(savedir,'Apr02_sumCells_fixedTr'),'AllData','-v7.3')


% Test on new data
testFeatures

% Pass CNN image features to trained classifier
predictedLabels = predict(classifier, testFeatures, 'ObservationsIn', 'columns');

% Get the known labels
testLabels = testSet.Labels;




end
