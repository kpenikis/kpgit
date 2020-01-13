function savePopMPHdata
%
% plotNormMPHs
%
%  Original plots from PopulationVariability, now separate, with options
%  for sorting the units.
%
%  Intended to help categorize response types.
%
%  NOW only used to get and save MPH data, after changes to UnitData.
%

global AMrates binsmth

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
% spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

binsmth   = 20;
trMax     = 40;
getMPH    = 1;
Stimuli   = [];
Duration  = 500;
skipOnset = 0;


%%
%currently: common spike shift

FR_vec   = nan(numel(UnitData),500,5);
zFR_vec  = nan(numel(UnitData),500,5);
FR_Warn  = nan(numel(UnitData),500);
zFR_Warn = nan(numel(UnitData),500);
sp_trs   = nan(numel(UnitData),500,5,10);
% sp_Warn_trs  = nan(numel(UnitData),500,10);

for iUn = 1:numel(UnitData)
    
    %     spkshift = [];
    
    % - - - -   Get raw spike data   - - - - -
    get_trial_data_posthoc
    
    %     if isempty(spkshift)
    %         continue
    %     end
    
    % - - - -   Pull relevant data   - - - - -
    
    for stid = 2:6
        
        if isempty(MPH(MPH.ThisStimID==stid,:))
            continue
        end
        
        % Save 10 random trials
        thisRaster = vertcat(MPH(MPH.ThisStimID==stid,:).raster{:});
        ridx = randperm(size(thisRaster,1),10);
        sp_trs(iUn,1:size(thisRaster,2),stid-1,:) = thisRaster(ridx,:)';
        
        % Normalized FR within this unit, referenced to silence
        zFR_vec(iUn,1:size(thisRaster,2),stid-1) = mean(vertcat(MPH(MPH.ThisStimID==stid,:).zFR{:}),1);
        
        % Get smoothed FR
%         FR_vec(iUn,1:size(thisRaster,2),stid-1)  = mean(vertcat(MPH(MPH.ThisStimID==stid & MPH.PrevPd==AMrates(stid-1),:).FRsmooth{:}),1);
        FR_vec(iUn,1:size(thisRaster,2),stid-1)  = mean(vertcat(MPH(MPH.ThisStimID==stid,:).FRsmooth{:}),1);
        
    end 
    
    % Save avg Warn response for this unit
    FR_Warn(iUn,:)  = mean(Warn_FR,1);
    zFR_Warn(iUn,:) = mean(Warn_zFR,1);
    
end %iUn

savedir = fullfile(fn.figs,'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% Save MPH data
save(fullfile(savedir,'MPHdata_gauss'),'zFR_vec','zFR_Warn','FR_vec','FR_Warn','sp_trs','-v7.3')

clear binsmth
end