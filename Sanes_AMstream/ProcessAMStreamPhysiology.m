
%% Process raw data (creates Info, Phys, SoundData files)

pp_make_InfoPhys( BLOCK, subject, session ) 


%% Spike sort with UMS (creates Spikes file)

Spikes = pp_sort_session( subject, session );


%% Manual spike sorting

[subject,session,channel,Spikes] = pp_launch_manual_sort( subject, session , channel );

Spikes = pp_save_manual_sort( subject, session, Spikes, channel, spikes );


%% Rasters and initial analysis

ap_plot_rasters(subject, session, [channels], [clus] )


