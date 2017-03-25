function remove_field(subject,session)

% Load Data structure
fn = set_paths_directories(subject,session);
filename = sprintf('%s_sess-%s_Data',subject,session);
load(fullfile(fn.processed,subject,filename));

% Go through each cluster and call analysis programs

allclusters = fieldnames(Data);
for unit = allclusters'
    
    chname = sprintf( '%s_sess-%s_%s',subject,session,unit{:});
    load(fullfile(fn.sess_data,chname))
    eval(sprintf('%s = rmfield(%s,''blk'');',unit{:},unit{:}));
    save(fullfile(fn.sess_data,chname),unit{:},'-v7.3');
    
end