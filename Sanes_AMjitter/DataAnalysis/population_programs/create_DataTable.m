function create_DataTable(subject,METRIC)

close all
% Set up useful variables
DataTable = table;
%'Session' 'cluname' 'SU/MU' 'HP LP' 'dBSPL' 'CenterRate' 'BehavState' 'Depth' 'Jitter' 'TuningMean' 'dprime' 'indVar'
DataTable.Session    = 'placeholder';
DataTable.cluname    = 'placeholder';
DataTable.UnitType   = discretize(4,[0.5 1.5 2.5 3.5 4.5],'categorical',{'unk', 'SU', 'MU', 'noise'});
DataTable.Noiseband  = [0 0];
DataTable.dBSPL      = 0;
DataTable.CenterRate = 0;
DataTable.BehavState = categorical({'P'},{'P' 'D' 'A'},{'Passive' 'Drinking' 'Active Behavior'});
DataTable.depth      = 0;
DataTable.jitter     = 0;
DataTable.TuningMean = 0;
DataTable.dprimeFR   = 0;
DataTable.dprime25ms = 0;
DataTable.DiscrmVar  = 'placeholder';



% Get Data files
fn = set_paths_directories;
searchname = sprintf('%s_sess-*_Data.mat',subject);
alldatafiles = dir(fullfile(fn.processed,subject,searchname));
try
    for ii = 1:numel(alldatafiles)
        % Load Data file
        clear Data
        load(fullfile(fn.processed,subject,alldatafiles(ii).name));
        
        [~,session] = strtok(alldatafiles(ii).name,'-');
        session = session(2:3);
        fn = set_paths_directories(subject,session);
        
        % Load each clu file
        allclusters = fieldnames(Data);
        allclusters = allclusters(strncmp('ch',allclusters,2));
        
        for unit = allclusters'
            
            % Load cluster data file
            clear clustruct Data_clu
            cluefilename = sprintf( '%s_sess-%s_%s',subject,session,unit{:});
            clustruct = load(fullfile(fn.sess_data,cluefilename));
            clustruct = clustruct.(unit{:});
            
            for ibk = 1:numel(clustruct.block)
                for ip = 1:numel(clustruct.block(ibk).pars)
                    stimpars = clustruct.block(ibk).pars(ip).pars;
                    
                    % Get independent variable discriminated
                    if isfield(Data,'depth') && any(sum(Data.depth==[ibk ip]',1)==2)
                        indVar = 'depth';
                    elseif isfield(Data,'jitter') && any(sum(Data.jitter==[ibk ip]',1)==2)
                        indVar = 'jitter';
                    end
                    nmfieldname = ['nm_' indVar];
                    
                    
                    for id = 1:numel(clustruct.block(ibk).pars(ip).stimvals)
                        
                        this_depth = clustruct.block(ibk).pars(ip).stimvals(id).depth;
                        these_jitters = clustruct.block(ibk).pars(ip).stimvals(id).jitter;
                        
                        %%% move behav here
                        behavs = cellfun(@unique,clustruct.block(ibk).pars(ip).stimvals(id).behav,'UniformOutput',false);
                        for ib = 1:numel(behavs)  %separate behavioral states
                            
                            match_behav = strcmp(behavs{ib},clustruct.block(ibk).pars(ip).stimvals(id).behav);
                            these_TuningMeans = clustruct.block(ibk).pars(ip).(METRIC)(id).mean(match_behav);
                            
                            for ij = 1:numel(these_jitters)
                                
                                this_jitter = char(these_jitters(ij));
                                this_TuningMean = these_TuningMeans(ij);
                                
                                % neurometric data is indexed in a different
                                % way
                                if ~isempty(clustruct.block(ibk).pars(ip).(nmfieldname))
                                    switch indVar
                                        case 'depth'
                                            condVar  = this_jitter;
                                            discrVar = [num2str(this_depth) '_depth'];
                                        case 'jitter'
                                            condVar  = [num2str(this_depth) '_depth'];
                                            discrVar = this_jitter;
                                    end
                                    
                                    if (strcmp(indVar,'jitter') && strcmp(this_jitter,'0')) || (strcmp(indVar,'depth') && strcmp(discrVar,'0_depth'))
                                        %if this is the NOGO value
                                        this_dprime = nan;
                                        
                                    else
                                        idxCls = find([clustruct.block(ibk).pars(ip).(nmfieldname).classifier.binsize]==25);
                                        idxOut = find(strcmp(behavs{ib},{clustruct.block(ibk).pars(ip).(nmfieldname).classifier(1).output.behav}));
                                        idx1 = find(strcmp({clustruct.block(ibk).pars(ip).(nmfieldname).stim{:,1}},condVar));
                                        idx2 = find(strcmp(clustruct.block(ibk).pars(ip).(nmfieldname).stim{idx1,2},discrVar));
                                        
                                        this_dprime = clustruct.block(ibk).pars(ip).(nmfieldname).classifier(1).output(idxOut).dprime.data{idx1}(idx2-1,2);
                                        this_dprime_25ms = clustruct.block(ibk).pars(ip).(nmfieldname).classifier(idxCls).output(idxOut).dprime.data{idx1}(idx2-1,2);
                                    end
                                    
                                else
                                    this_dprime = nan;
                                end
                                
                                
                                % Add new row to DataTable
                                %'Session' 'cluname' 'SU/MU' 'Noise band' 'dBSPL' 'CenterRate' 'BehavState' 'Depth' 'Jitter' 'TuningMean' 'dprime' 'indVar'
                                DT_addrow = {session unit{:} ...
                                    discretize(clustruct.labels(1,2),[0.5 1.5 2.5 3.5 4.5],'categorical',{'unk', 'SU', 'MU', 'noise'})...
                                    [stimpars(1) stimpars(2)] stimpars(3) stimpars(4) ...
                                    categorical(behavs(ib),{'P' 'D' 'A'},{'Passive' 'Drinking' 'Active Behavior'})...
                                    this_depth this_jitter this_TuningMean this_dprime this_dprime_25ms indVar};
                                
                                DataTable = [DataTable; DT_addrow];
                                
                                
                            end
                        end
                    end
                    
                end
            end
            
        end
        
    end
    
catch
    keyboard
end

keyboard


%Export dataset
ds = dataset(daymat,behavior,neural);
export(ds,'file',[statsdirectory,fname,'_thresholds.csv'],'delimiter','comma');

% for each subplot
% get relevant datapoints and plot



end


