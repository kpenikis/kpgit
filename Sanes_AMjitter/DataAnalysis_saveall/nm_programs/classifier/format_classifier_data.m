function [this_output,stim,cludata,Data_clu] = format_classifier_data(cludata,is,ip,METRIC,binsize,indVar,Data_clu)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each
%

global fn subject session channel paramdir

% Get struct of this stimulus
stimdata = cludata.block(is).pars(ip);

% Get raster indices for this stimulus
if ~isfield(stimdata.stimvals,'raster_idxs')
% if strcmp('FR',METRIC)
    fprintf('   (getting raster idxs)\n')
    [stimdata,raster_out] = get_raster_indices(stimdata,Data_clu.raster,cludata.block(is).block,indVar);
    cludata.block(is).pars(ip) = stimdata;
    Data_clu.raster = raster_out;
end


%% Run classifier only on identical behavioral state

this_output=struct();

behavstates = unique([stimdata.stimvals.behav]);

for ibs = 1:numel(behavstates)
    
    bs = behavstates(ibs);
    
    % Make string for title and savename
    [title_str,savename] = make_title_savename_str(stimdata, cludata.block(is).block, bs,...
        subject,session,channel,cludata.labels(1,1),METRIC,binsize,indVar,'classif');
    
    % Get raster indices that correspond to this behavioral state
    ridx_bs = [stimdata.stimvals.raster_idxs];
    ridx_bs = ridx_bs( strcmp(bs{:},[stimdata.stimvals.behav]) );
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Get classifier data
    output.behav = bs{:};
    try
    [PYdata,dprime,stim] = get_classifier_data( Data_clu.raster(ridx_bs), METRIC, binsize, indVar);
    catch
        aaa=234;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output.PYdata.data = PYdata;
    output.dprime.data = dprime;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Fit classifier data
    [fitdata,hF] = fit_classifier_data( output, stim, title_str, indVar );
    output.PYdata.fitdata = fitdata.PYdata;
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

