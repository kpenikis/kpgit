function updateClassifierResStimData
% 

global AMrates rateVec_AC rateVec_DB trMin 

%% 
% Load Unit files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------
[sigUnits,UnitData] = identifyResponsiveUnits(UnitData);

% Load classifier results
load(fullfile(fn.processed,'MPHclassifier','ClassData'));
q=load(fullfile(fn.processed,'MPHclassifier','ClassData_shuff'));
DataSh = q.Data;
clear q

% Other settings
histbinsize = 0.025;
trMin       =  10;

AMrates     = [2 4 8 16 32];

q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;


%%  Call <<getFRhist>>

for iUn = 1:size(Data,1)
    %-----------------------------------
    % Get trial-by-trial firing history
    Data(iUn,:)   = getFRhist( Data(iUn,:),   UnitData(iUn), spkshift );
    
    %DataSh(iUn,:) = getFRhist( DataSh(iUn,:), UnitData(iUn), spkshift );
    % can't run for shuffled data! because classifier shuffled spike times
    % are not saved. instead, running classifier again.
    
    %-----------------------------------
end %iUn


% Resave Data structures
savedir = fullfile(fn.processed,'MPHclassifier');
save(fullfile(savedir,'ClassData'),'Data','-v7.3')



end


