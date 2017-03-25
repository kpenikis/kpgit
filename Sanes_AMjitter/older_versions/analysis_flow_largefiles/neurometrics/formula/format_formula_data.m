function [this_output,cludata] = format_formula_data(cludata,is,ip,indVar)
% 

global subject session channel xlimits

% Get struct of this stimulus
stimdata = cludata.analyses.stim(is).pars(ip);

% Get raster indices for this stimulus
if ~isfield(stimdata.depths(1),'raster_idxs')
    fprintf('getting raster idxs')
    stimdata = get_raster_indices(stimdata,cludata);
    cludata.analyses.stim(is).pars(ip) = stimdata;
end


%% Separate by behavioral state

this_output=struct();

behavstates = unique([stimdata.depths.behav]);

for ibs = 1:numel(behavstates)
    
    bs = behavstates(ibs);
    
    % Make string for title and savename
    [title_str,savename] = make_title_savename_str(cludata.analyses.stim(is).pars(ip),bs,...
        subject,session,channel,cludata.label(1),indVar);
    
    % Get raster indices that correspond to this behavioral state
    ridx_bs = [stimdata.depths.raster_idxs];
    ridx_bs = ridx_bs(strcmp([stimdata.depths.behav], bs ));
    
    % Get minimum stimulus duration
    dur = min([cludata.raster(ridx_bs).stimDur]);
    dur = 'fulldur';

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Create FR matrix
    [FR_mat,stim,ntrials] = format_FRmat( cludata.raster(ridx_bs), dur, indVar );
    
    % Calculate dprime based on the formula 
    dprime_mat = calculate_dprime_formula( FR_mat );
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Format output for cludata struct
    output.behav = bs;
    output.stim = stim;
    output.FR_mat = FR_mat;
    output.dprime_mat = dprime_mat;

    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Fit with sigmoid and calculate dprime threshold
    fitdata = fit_calc_dprime_threshold( dprime_mat, ntrials );
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output.fitdata = fitdata;
    
    this_output(ibs).output = output;
    
    
    %% Plot and save figure
    
    % Set up formatting
    [~, plotOptions] = setOptions;
    switch indVar
    case 'depth'
        colors = copper(numel(output.dprime_mat));
    case 'jitter'
        colors = winter(numel(output.dprime_mat));
        plotOptions.xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
    end
    legtext=cell(1,size(output.stim,1));
    for ic=1:size(output.stim,1)
        legtext{ic}=char(output.stim{ic,1});
    end
    
    % Set up figure
    hF = figure; hold on
    scrsz = get(0,'ScreenSize');
    set(hF,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    hS=subplot(1,1,1);
    
    % Plot
    hL = plot_neurometric( output.dprime_mat, fitdata, plotOptions, hS, xlimits, colors );
    
    % Finish formatting
    legend(hL,legtext,'Interpreter','none')
    title(hS,'dprime by FR formula fit with sigmoid (nlinfit)')
    suptitle(title_str)
    
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



function [title_str,savename] = make_title_savename_str(parstruct,bs,subject,session,channel,clu,indVar)

% Get pars for title and savename
block = parstruct.block;
HP    = parstruct.pars(1);
LP    = parstruct.pars(2);
dB    = parstruct.pars(3);
rate  = parstruct.pars(4);

% Add stimulus info in title
title_str = sprintf('ch %i clu%i:  neurometric %s - FR formula  |  %s  |  noise: %i - %i Hz  |  %idB  |  AM %i Hz  |  blk%i',...
    channel,clu,indVar,bs{:},HP,LP,dB,rate,block);

savename  = sprintf('%s_%s_nm%s_formulaFR_ch%i_clu%i_AM%iHz_%idB_%i-%i_blk%i',...
    subject,session,indVar,channel,clu,rate,dB,HP,LP,block);

end
    


