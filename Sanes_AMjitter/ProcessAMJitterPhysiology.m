

subject   =  'IIIf_230115' ;
session   =  'LA' ;

channel   =  1 ;
clu       =  1 ;

fn = set_paths_directories(subject,session);


%% Process raw data (creates Info, Stim, and Phys files)

pp_prepare_format( BLOCKS, subject, session ) 


%% Spike sort with UMS (creates Spikes file)

Spikes = pp_sort_session( subject, session );


%% Manual spike sorting

[subject,session,channel,Spikes] = pp_launch_manual_sort( subject, session , channel );

Spikes = pp_save_manual_sort( subject, session, Spikes, channel, spikes );


%% Create data structure, raster/psth plots (creates Data file)

create_Data_struct( subject, session )

    %raster = pp_plot_rasters_wStim(subject,session,channel,clu);

    
    
%% Run analysis plots

anResponsePopulationPlots( subject )


    
%% Neurometric analyses

Data = getALL_nmData(subject,session);

    %cludata = ap_nmData_depth(cludata)

    
    
%% Collect Standard Period Responses from the population

StandardPd_Responses( subject )


% Plot various measures of activity during standard period

standardPd_anyPlot( subject, plottype, var1, var2, categorization)

% ('IIIf_230115', psth,   '',   '',   'ALL')
% ('IIIf_230115', psth,   '',   '',   'Corr')
% ('IIIf_230115', psth,   '',   '',   'VS')
% ('IIIf_230115', psth_diff,  '',   '',   'ALL')
% ('IIIf_230115', comparison, 'FR',         '', 'ALL')
% ('IIIf_230115', comparison, 'FRnorm',     '', 'ALL')
% ('IIIf_230115', comparison, 'PeakTime',   '', 'ALL')
% ('IIIf_230115', comparison, 'PeakTrough', '', 'ALL')
% ('IIIf_230115', regression, 'FRidx', 'FR100ms_idx', 'ALL')
% ('IIIf_230115', regression, 'FRidx', 'FR250ms_idx', 'ALL')
% ('IIIf_230115', regression, 'FRidx', 'prevRate',    'ALL')


%%  - - - - OR - - - - - 

%% Plot SU psths and within unit metrics

ScreenUnitResponses( subject )







