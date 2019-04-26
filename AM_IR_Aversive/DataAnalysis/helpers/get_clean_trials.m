function [theseTDidx,Ntrials,minDur,allStim] = get_clean_trials(TrialData,artifact_trs,spl,lpn)
% called by ap_plot_rasters (aversive experiment)
% 


% Find trials that match the current spound parameters (spl, lpn)
idx_params = find(TrialData.SPL==spl & TrialData.LP==lpn);


if ismember('ITIflag',TrialData.Properties.VariableNames)
    
    % Find SAFE trials when animal was on spout for at least 90% of the time,
    % and all WARN trials
    idx_spout = find( (TrialData.ITIflag==0 & TrialData.trID>1 & TrialData.Spout>0.9) ...
        |  TrialData.trID==1 ...
        | (TrialData.ITIflag==1 & TrialData.Spout>0.75));
else
    idx_spout = 1:size(TrialData,1)';
end

% Find trials that are not flagged for artifact (not including Warn trials)
[~,~,rem] = intersect(find(TrialData.trID==1),artifact_trs);
artifact_trs(rem) = [];
idx_clean = 1:size(TrialData,1);
idx_clean ( artifact_trs ) = [];

% Get the indices of TrialData that comply with all of these filters
theseTDidx = intersect(intersect(idx_params,idx_spout),idx_clean');

TD = TrialData(theseTDidx,:);

% Get number of trials and duration of each stimulus
theseStim = unique(TD.trID);
Ntrials = nan(1,length(theseStim));
StimDur = nan(1,length(theseStim));
for id = 1:numel(theseStim)
    Ntrials(id) = sum(TD.trID==theseStim(id));
    StimDur(id) = mode(TD.offset(TD.trID==theseStim(id)) - TD.onset(TD.trID==theseStim(id)));
end
minDur = min(StimDur);


% And find what the stimuli are
allStim = unique(TrialData.trID(theseTDidx));

end

