
SplitUnits = struct;

SplitUnits(1).un1 = 'WWWf_SB_3_17';
SplitUnits(1).un2 = 'WWWf_SC_3_45';

SplitUnits(2).un1 = 'WWWf_YA_11_11';
SplitUnits(2).un2 = 'WWWf_YB_11_5';

fn = set_paths_directories('','',1);
save(fullfile(fn.processed,'SplitUnits'),'SplitUnits','-v7.3')

