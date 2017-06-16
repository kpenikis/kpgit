function add_an_datapoints_table( session, cluname, clulabel,...
                        stimpars, behav, indVar, dp)

global DataTable

try
for istim = 1:numel(dp.depth)
    
    this_jitter = char(dp.jitter(istim)); %char(stimvals(istim,1));
    this_depth  = dp.depth(istim);
    
    %'Session' 'cluname' 'SU/MU' 'HP LP' 'dBSPL' 'CenterRate' 'BehavState' 'Depth' 'Jitter' 'dprime'... 'indVar'
    DT_addrow = {session cluname ...
        discretize(clulabel, [0.5 1.5 2.5 3.5 4.5], 'categorical',{'unk', 'SU', 'MU', 'noise'})...
        stimpars(1) stimpars(2) stimpars(3) stimpars(4) ...
        categorical(behav, {'P' 'D' 'A'}, {'Passive' 'Drinking' 'Active Behavior'})...
        this_depth this_jitter dp.baselineFR(istim) indVar };
    
    dpfns = fieldnames(dp);
    
    for ifn = dpfns'
        
        if any(strcmp(ifn{:},{'depth' 'jitter' 'baselineFR'}))
            continue
        end
                
        DT_addrow{end+1}  = dp.(ifn{:})(istim);
        
    end
    
    if size(DataTable,2) ~= size(DT_addrow,2)
        keyboard
    end
    
    try
    DataTable = [DataTable; DT_addrow];
    catch
        keyboard
    end
    
end
catch
    keyboard
end
end



