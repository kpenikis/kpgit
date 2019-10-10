
%% Plot MPHs of Population

plotPopNormMPH  %optionally creates mph datafile
    get_trial_data_posthoc
    
clusterMPHtypes %just plots
    mph_pca


quickPlotPopNormSpeech  % no mph, so onsets of stim
    get_trial_data_speech



%% Compare activity in matched MPs

GetMatchedMPData
    makeMPHtable


MPcontext_wins %compare context discrimination results from full period and 62 ms


%% Random stuff 

NunitsPerSess  %histogram of # units per session

