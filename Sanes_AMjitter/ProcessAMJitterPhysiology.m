

subject   =  'IIIf230115' ;
session   =  'LA' ;

channel   =  1 ;
clu       =  1 ;


%% Process raw data (creates Info, Stim, and Phys files)

pp_prepare_format( BLOCKS, subject, session ) 


%% Spike sort with UMS (creates Spikes file)

Spikes = pp_sort_session( subject, session );


%% Manual spike sorting

[subject,session,channel,Spikes] = pp_launch_manual_sort( subject, session , channel );

Spikes = pp_save_manual_sort( subject, session, Spikes, channel, spikes );


%% Create data structure, raster/psth plots (creates Data file)

create_Data_struct( subject, session )


%% Run analysis plots

runALL_anPlots( subject, session );

