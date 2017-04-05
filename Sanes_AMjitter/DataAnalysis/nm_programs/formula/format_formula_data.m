function [this_output,stim,raster] = format_formula_data(raster,indVar)
% 

global subject session channel clu xlimits fn paramdir


%% Separate by behavioral state

this_output=struct();

behavstates = unique({raster.behaving});

for ibs = 1:numel(behavstates)
    
    bs = behavstates{ibs};
    
    % Make string for title and savename
    [title_str,savename] = make_title_savename_str(raster, bs,...
        subject,session,channel,clu,'',[],indVar,'formulaFR');
    
    % Get raster indices that correspond to this behavioral state
    ridx_bs = find( strcmp({raster.behaving}, bs) );
    
    
    % Get minimum stimulus duration
%     dur = min([Data_clu.raster(ridx_bs).stimDur]);
    dur = 'fulldur';

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Create FR matrix
    [FR_mat,stim,ntrials] = format_FRmat( raster(ridx_bs), dur, indVar );
    
    % Calculate dprime based on the formula 
    dprime_mat = calculate_dprime_formula( FR_mat );
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Format output for cludata struct
    output.behav = bs;
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
            ALLcolors = winter(numel(plotOptions.colSelect));
            for ic=1:size(stim,1)
                colors(ic,:) = ALLcolors(strcmpi(stim{ic,1},plotOptions.colSelect),:);
            end
            plotOptions.xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
    end
    legtext=cell(1,size(stim,1));
    for ic=1:size(stim,1)
        legtext{ic}=char(stim{ic,1});
    end
    
    % Set up figure
    hF = figure; hold on
    scrsz = get(0,'ScreenSize');
    set(hF,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    hS=subplot(1,1,1);
    
    % Plot
    hL = plot_neurometric( output.dprime_mat, fitdata, plotOptions, hS, xlimits, colors );
    
    % Finish formatting
    legend(hL,legtext,'Interpreter','none','Location','best')
    title(hS,'dprime by FR formula fit with sigmoid (nlinfit)')
    suptitle(title_str)
    
    
    %% Finish and save figure

    %%%%%%%%%%%%
    %   Save   %
    %%%%%%%%%%%%
    
    % Set save directory
    sv_dir = fullfile(fn.nmplots,indVar,'formula',paramdir);
    if ~exist(sv_dir,'dir')
        mkdir(sv_dir)
    end
    
    % Save figure
    set(hF,'PaperOrientation','landscape');
    print(hF,'-dpdf',fullfile(sv_dir,savename),'-bestfit')
    
    
    
    
end



end




