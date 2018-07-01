

analyzeMPH('FR','iBMF_FR')

% after removing non-sync points                                            *** talk to dan about significance
analyzeMPH('VS','iBMF_FR')

% after adding correlation to stats
% also plot normalized context comparisons 
% (for rate/dist BMF, also prev pd and prior FR might as well)
analyzeMPH('FR','all')
analyzeMPH('VS','all')  %                                                   *** FIND OUT IF VS VALUES DIFFERENT IN CONTEXT PROGRAM


% after saving avg IR obs plots
predObs_population('FR')

% after adding VS prediction (what to do with non-sync rates?)              *** too complicated, even in calculating observed response... skip
predObs_population('VS')


% INTERMEDIATE ANALYSES



% Write code to analyze power of spiketimes while drinking vs not



