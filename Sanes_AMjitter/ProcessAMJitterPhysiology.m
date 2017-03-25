

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

getALL_anData( subject, session );

    %ap_xJitter(subject,session,channel,clu,METRIC)
    
    
%% Neurometric analyses

Data = getALL_nmData(subject,session);

    %cludata = ap_nmData_depth(cludata)


