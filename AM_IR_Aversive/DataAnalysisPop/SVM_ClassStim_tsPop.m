function SVM_ClassStim_tsPop
% SVM_ClassStim_tsPop
%
%  SVM classifier for segments of Pdc and Irr stimuli.
%
%
%  KP, 2019-12
%


% close all

whichIrr    = 'AC';
StimDur     = 500;
convwin     = 10;
convwin     = ones(1,convwin).*(1/convwin);


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');

% Load spikes data (created in cumulativeSpikeCount)
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

% Find just Apr02 cells
sess_idx = find(strcmp(UnitInfo.Subject,'AAB_265054') & strcmp(UnitInfo.Session,'Apr02-AM'));
Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(intersect(Un_Indices,sess_idx),:,:,:);


% rows: observations
% columns: dimensions
% table: trial * stim  X  cell 
%   time?

%% Data: avg over time

TrainingData = [];
CorrResponse = [];

for ist = 1:size(Cell_Time_Trial_Stim,4)
    
    nTr = find(~isnan(mean(mean(Cell_Time_Trial_Stim(:,:,:,ist),2),1)),1,'last');
    
    % automatically skips empty stim
    for itr = 1:nTr
        TrainingData = [ TrainingData; sum(Cell_Time_Trial_Stim(:,:,itr,ist),2)' ];
        CorrResponse = [ CorrResponse; ist ];
    end
    
end


AllData = [TrainingData CorrResponse];

savedir = fullfile(fn.figs,'SVMclass','testing');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
save(fullfile(savedir,'Apr02_sumTime'),'AllData','-v7.3')




% Classify using avg N spks across cells
% Linear SVM:
%ovo   53.8%
%ova   50.4%
% Basic linear classifier:
%      48.4%

% Classify using timeseries avg across cells
% Linear SVM (1 ms):
%ovo   43.6%
%ova   44.2%
% Basic linear classifier:
%      15.6% (1 ms)
%      40.6% (10ms smoothed)



%% Data: avg over cells

TrainingData = [];
CorrResponse = [];
AllData      = [];

for ist = 1:size(Cell_Time_Trial_Stim,4)
    
    nTr = find(~isnan(mean(mean(Cell_Time_Trial_Stim(:,:,:,ist),2),1)),1,'last');
    
    % automatically skips empty stim
    for itr = 1:nTr
        TrainingData = [ TrainingData; mean(Cell_Time_Trial_Stim(:,:,itr,ist),1) ];
        CorrResponse = [ CorrResponse; ist ];
    end
    
end

AllData = [TrainingData CorrResponse];

savedir = fullfile(fn.figs,'SVMclass','testing');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
save(fullfile(savedir,'Apr02_sumCells'),'AllData','-v7.3')

end