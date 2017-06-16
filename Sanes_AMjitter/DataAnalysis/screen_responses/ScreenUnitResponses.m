function ScreenUnitResponses(subj)

global fn subject only75
subject=subj;

fn = set_paths_directories;

%*************
only75 = 1;
%*************

close all


%% FOR EACH INDEPENDENT VARIABLE

for iV = { 'jitter' 'depth'  }

    indVar = iV{1};
    nmfieldname = ['nm_' indVar];

%% Get all processed Data files
searchname = sprintf('%s_sess-*_Data.mat',subject);
alldatafiles = dir(fullfile(fn.processed,subject,searchname));

for ii = 1:numel(alldatafiles)
    
    % Load Data file
    clear Data
    load(fullfile(fn.processed,subject,alldatafiles(ii).name));
    
    % Get this session label and set directories
    [~,session] = strtok(alldatafiles(ii).name,'-');
    session = session(2:3);
    fn = set_paths_directories(subject,session);
    
    % Load each cluster file
    allclusters = fieldnames(Data);
    allclusters = allclusters(strncmp('ch',allclusters,2));
    
    for unit = allclusters'
        
        % *****  SKIP IF << create_Data_struct >> (new version) NOT RUN YET  ***** %
        if ~isfield(Data.(unit{:}),'stimdata'), continue, end
        % ************************************************************************ %
        
        % ONLY SU FOR NOW
        if Data.(unit{:}).labels(2) ~= 2
            continue
        end
        
        % Check if this session contains the discriminated variable
        stimdata = Data.(unit{:}).stimdata;
        get_sess_data = [];
        for is = 1:numel(stimdata)
            if isfield(stimdata(is).pars,nmfieldname)
                get_sess_data = [get_sess_data is];
            end
        end
        
        if isempty(get_sess_data), continue, end
        
        %% If there is relevant data in this session, load cluster data file
        
        % Otherwise, get that data
        for is = get_sess_data
            for ip = 1:numel(stimdata(is).pars)                
                
                % Get stimulus parameters and behavioral state
                stimpars = stimdata(is).pars(ip).pars;
                behavs = unique(stimdata(is).pars(ip).stimvals(:,3));
                
                for ib = 1:numel(behavs)  %separate behavioral states
                    
                    %% Get these datapoints
                    
                    rasters = stimdata(is).pars(ip).(nmfieldname);
                    rasters = rasters(strcmp({rasters.behaving},behavs(ib)));
                    
                    %% Plot PSTH of each jitter with periodic
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    standardPd_PSTHcomparisons( rasters,...
                        session, unit{:}, Data.(unit{:}).labels(1,2),...
                        stimpars, behavs(ib), indVar, is )
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    close all
                    
                end %behavs
            end %pars
        end %stim block
    end %unit
end %session
end % for indVar


end





