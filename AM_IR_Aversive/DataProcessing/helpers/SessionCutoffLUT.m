fn = set_paths_directories;

% 2019-08-27
T=table;

T.subject   = 'AAB_265054';
T.session   = 'Mar25-AM';
T.ch        = nan;
T.clu       = 598;
T.CutTrial  = nan;

T(end+1,:) = {'AAB_265054' 'Mar25-AM' nan 632 nan};
T(end+1,:) = {'AAB_265054' 'Mar26-AM' nan 1104 nan};
T(end+1,:) = {'AAB_265054' 'Mar26-AM' nan 1137 nan};
T(end+1,:) = {'AAB_265054' 'Mar26-AM' nan 1463 nan};
T(end+1,:) = {'AAB_265054' 'Mar26-AM' nan 1582 nan};
T(end+1,:) = {'AAB_265054' 'Mar28-AM' nan 455 nan};


% Add to table
load(fullfile(fn.processed,'UnLUTcutDB.mat'),'-mat')

T(end+1,:) = {'AAB_265054' 'Mar28-AM' nan 455 nan};


% Save
save(fullfile(fn.processed,'UnLUTcutDB'),'T','-v7.3')

