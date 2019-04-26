
%%  KS PROCESSING

ProcessPhys_SynKS( 'Aug25-AM' )

% MANUAL SORTING 
%    cd ../sorting
%    source activate phy
%    phy template-gui params.py

KS_to_Spikes  % <-- MUST run this after manual sorting


%% Rasters and initial analysis 

% To check single session:
allRastersSession(SUBJECT,SESSION)
screenUnitsRasters(SUBJECT,SESSION)
AssessUnits(SUBJECT,SESSION) %be sure to change save name in file
getWaveformShapes('Units_filename','Plot_savename')

% Add to overall dataset: 
AssessUnits
assessWaveformShapes('Units','Plot_savename') 
    %merged units weird now
    %adds to UnitInfo

predObserved  
    %not set up for merged units
pdcRespRepetition
MPHcontextComparisons



%%  UMS PREPROCESSING  *use _UMS dir

% Set input vars
subject = 'WWWf_253400';
session = 'PA';
BLOCK = 41;

% Process raw data (creates Info, Phys, TrialData files)
ProcessPhys_UMS( BLOCK, subject, session ) 

% Spike sort with UMS (creates Spikes file)
Spikes = pp_sort_session( subject, session );
keyboard

% Manual spike sorting
[subject,session,channel,Spikes] = pp_launch_manual_sort( subject, session , channel );
Spikes = pp_save_manual_sort( subject, session, Spikes, channel, spikes );


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


%% Rasters and initial analysis (old)

AssessUnits
edit SplitUnits.mat
plotUnitDistributions

plotRasters_allUnits
plotRasters_eachUnit

predObserved

DirectPdComparison
    makeMPHtable
    plotMPHs


% Other analyses
MPHpolarHistograms
OnsetWinRatio
FullCycleRatio
FullCycleFF




%%  ARCHIVE: Rasters and initial analysis

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

