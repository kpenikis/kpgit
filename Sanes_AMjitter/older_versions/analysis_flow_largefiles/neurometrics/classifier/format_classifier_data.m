function [this_output,cludata] = format_classifier_data(cludata,is,ip,METRIC,binsize,indVar)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each
%

global subject session channel

% Get struct of this stimulus
stimdata = cludata.analyses.stim(is).pars(ip);

% Get raster indices for this stimulus
if ~isfield(stimdata.depths(1),'raster_idxs')
    fprintf('   (getting raster idxs)\n')
    stimdata = get_raster_indices(stimdata,cludata);
    cludata.analyses.stim(is).pars(ip) = stimdata;
end


%% Run classifier only on identical behavioral state

this_output=struct();

behavstates = unique([stimdata.depths.behav]);

for ibs = 1:numel(behavstates)
    
    bs = behavstates(ibs);
    
    % Make string for title and savename
    [title_str,savename] = make_title_savename_str(cludata.analyses.stim(is).pars(ip),bs,...
        subject,session,channel,cludata.label(1),METRIC,binsize,indVar);
    
    % Get raster indices that correspond to this behavioral state
    ridx_bs = [stimdata.depths.raster_idxs];
    ridx_bs = ridx_bs(strcmp([stimdata.depths.behav], bs ));
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Get classifier data
    output.behav = bs;
    output = get_classifier_data( cludata.raster(ridx_bs), METRIC, binsize, indVar);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Fit classifier data
    [fitdata,hF] = fit_classifier_data( output, title_str, indVar );
    output.fitdata = fitdata;
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Insert classifier output and raster indices into cludata struct
    this_output(ibs).output = output;
    
    
    %% Finish and save figure
            
    % Set save directory
    datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
    sv_dir = fullfile(datadir,subject,'^nm_plots',session,indVar);
    if ~exist(sv_dir,'dir')
        mkdir(sv_dir)
    end
    
    % Save figure
    set(hF,'PaperOrientation','landscape');
    print(hF,'-dpdf',fullfile(sv_dir,savename),'-bestfit')
    
    
    
    
end


end

function [title_str,savename] = make_title_savename_str(parstruct,bs,subject,session,channel,clu,METRIC,binsize,indVar)

% Get pars for title and savename
block = parstruct.block;
HP    = parstruct.pars(1);
LP    = parstruct.pars(2);
dB    = parstruct.pars(3);
rate  = parstruct.pars(4);

% Add stimulus info in title
title_str = sprintf('ch %i clu%i:  neurometric %s - %s bin %i  |  %s  |  noise: %i - %i Hz  |  %idB  |  AM %i Hz  |  blk%i',...
    channel,clu,indVar,METRIC,binsize,bs{:},HP,LP,dB,rate,block);

savename  = sprintf('%s_%s_nm%s_classif_%s-bin%i_ch%i_clu%i_AM%iHz_%idB_%i-%i_blk%i',...
    subject,session,indVar,METRIC,binsize,channel,clu,rate,dB,HP,LP,block);

end
    


