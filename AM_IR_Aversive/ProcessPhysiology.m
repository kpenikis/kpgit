
%%  AUTOMATED PREPROCESSING

% Set input vars
subject = 'WWWlf_253395';
% subject = 'WWWf_253400';
session = 'PA';
BLOCK = 41;

% Process raw data (creates Info, Phys, TrialData files)
pp_processPhys_UMS( BLOCK, subject, session ) 

% Spike sort with UMS (creates Spikes file)
Spikes = pp_sort_session( subject, session );


%% Re-make TrialData to include more ITI trials

SKIP_PHYS = 1;

% Set input vars
Subjects = {'WWWf_253400'};
Sessions = {'GA' 'KA' 'PA' 'QB' 'SA' 'SB' 'SC'};
BLOCKS   = [ 17   30   54   60   66   69   72 ];

for subj = Subjects
    for is = 1:numel(Sessions)
        
        subject = subj{:};
        session = Sessions{is};
        block = BLOCKS(is);
        
        % Process raw data (creates Info, Phys, TrialData files)
        pp_processPhys_UMS( block, subject, session, SKIP_PHYS )
        
        % Also redo artifact flagging, because trial numbers are different
        coldcall_identify_artifact( subject, session )
        
    end
end


%% Manual spike sorting

[subject,session,channel,Spikes] = pp_launch_manual_sort( subject, session , channel );

Spikes = pp_save_manual_sort( subject, session, Spikes, channel, spikes );


%% Rasters and initial analysis

ap_plot_rasters(subject, session, [channels], [clus] )

AssessResponses

predObs_population
Ztransitions
FFtransitions
cc_psth_rms
ContextMTFs
analyzeMPH



% ap_transitions
% ap_mph
% ap_zscore_plots
% ap_rcorr


%% Population level analyses
ap_zFR_population


%% ICAC Poster

ap_zFR_population
ap_mph_matched
ap_IRtrans_diffspks_v2

