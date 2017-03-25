function updatePopulationPlots(subject,indVar,METRIC)

close all
% Set up useful variables
figsStruct = struct;
DataTable = table;
% Session cluname HP LP dBSPL CenterRate BehavState Depth Jitter TuningMean dprime 


% Get Data files
fn = set_paths_directories;
searchname = sprintf('%s_sess-*_Data.mat',subject);
alldatafiles = dir(fullfile(fn.processed,subject,searchname));

for ii = 1:numel(alldatafiles)
    % Load Data file
    clear Data
    load(fullfile(fn.processed,subject,alldatafiles(ii).name));
    
    % Check for depth or jitter detection
    if ~isfield(Data,indVar)
        continue
    end
    
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
        channel = clustruct.labels(1,3);
        Data_clu = Data.(unit{:});
        
        iV_coords = Data.(indVar);
        
        for is=1:size(iV_coords,2)
            
            % Get string of stimulus parameters
            stimpars = clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).pars;
            behavs = cellfun(@unique,{clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals(end).behav});
            
            for ib = 1:numel(behavs)  %separate behavioral states
                
                parname = sprintf('pars%i_%i_%i_%i_%s',stimpars(1),stimpars(2),stimpars(3),stimpars(4),behavs{ib});
                
                if ~isfield(figsStruct,parname)
                    fighandle = ['hF' num2str(numel(fieldnames(figsStruct))+1)];
                    eval(sprintf('%s = figure;',fighandle))
                    figsStruct.(parname) = fighandle;
                end
                
                % Set current figure
                eval(sprintf('figure(%s)',figsStruct.(parname)))
                hold on
                
                %% Get data to plot
                switch indVar
                    case 'depth'
                        xvals  = [clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals.depth];
                        xvals = convert_depth_proptodB(xvals);
                        nconds = max(cellfun(@numel,{clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals.jitter}));
                        conds = unique(clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals(end).jitter,'stable');
                        
                        tuning_data = get_tuning_data(clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).(METRIC),nconds)
                        
                        
%                         clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).nm_depth.classifier(1).output.dprime.data{1};
%                         stimvals = clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals;
                        
                    case 'jitter'
                        % Get x axis values (jitter harder than depth)
                        %also need separate subplots
                        xvals = double(strtok(clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals(end).jitter,'_'));
                        xvals = convert_jitter_rangeHz(xvals,stimpars(4));
                        nconds = numel(clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals);
                        conds = cell(1);
                        for ic = 1:nconds
                            conds{ic} = [num2str(clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).stimvals(ic).depth) '_depth'];
                        end
                        
                        stimvals.(indVar)
                        
                        clustruct.block(iV_coords(1,4)).pars(iV_coords(2,4)).nm_jitter
                end
                
                tuningDATA = [clustruct.block(iV_coords(1,is)).pars(iV_coords(2,is)).(METRIC).mean];
                
            end
        end
        
        keyboard
        
    end
    
end

keyboard

%Export dataset
ds = dataset(daymat,behavior,neural);
export(ds,'file',[statsdirectory,fname,'_thresholds.csv'],'delimiter','comma');

% for each subplot
% get relevant datapoints and plot



end


