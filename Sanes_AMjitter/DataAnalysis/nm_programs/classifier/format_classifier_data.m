function [this_output,stim] = format_classifier_data(raster,METRIC,binsize,indVar)
% output = format_classifier_data(raster,METRIC,binsize,indVar)
%
%  First organizes all unique stimuli, then steps through to get AM depth
%  detection discriminability data for each
%

global fn subject session channel clu paramdir

this_output=struct();

%% Run classifier only on identical behavioral state
behavstates = unique({raster.behaving});

for ibs = 1:numel(behavstates)
    
    bs = behavstates{ibs};
    
    % Make string for title and savename
    [title_str,savename] = make_title_savename_str(raster, bs,...
        subject,session,channel,clu,METRIC,binsize,indVar,'classif');
    
    % Get raster indices that correspond to this behavioral state
    ridx_bs = find( strcmp({raster.behaving}, bs) );
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Get classifier data
    output.behav = bs;
    try
    [PYdata,dprime,stim] = get_classifier_data( raster(ridx_bs), METRIC, binsize, indVar);
    catch
        keyboard
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output.PYdata.data = PYdata;
    output.dprime.data = dprime;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Fit classifier data
    [fitdata,hF] = fit_classifier_data( output, stim, title_str, indVar );
%     output.PYdata.fitdata = fitdata.PYdata;
    output.dprime.fitdata = fitdata.dprime;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Insert classifier output and raster indices into cludata struct
    this_output(ibs).output = output;
    
    
    %% Finish and save figure
    
    %%%%%%%%%%%%
    %   Save   %
    %%%%%%%%%%%%
    
    % Set save directory
    sv_dir = fullfile(fn.nmplots,indVar,'classif',paramdir);
    if ~exist(sv_dir,'dir')
        mkdir(sv_dir)
    end
    
    % Save figure
    set(hF,'PaperOrientation','landscape');
    print(hF,'-dpdf',fullfile(sv_dir,savename),'-bestfit')
    
    
    
    
end


end

