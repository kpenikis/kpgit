function StandardPd_Responses(subj)

global InfoTable AmalgamatedSpikes OtherPdSpikes fn subject only75
subject=subj;

fn = set_paths_directories;

%*************
only75 = 1;
%*************

%%
if nargin<2 && ~exist('METRICS','var')
%****************************************%
METRICS = {'mean' 'std' 'p' 'FR' 'VS' 'meanphase'};
%****************************************%
end

%% Set up data and info

AmalgamatedSpikes = zeros(60000,250);
OtherPdSpikes = zeros(60000,250,5);

InfoTable = table;
%'Session' 'cluname' 'SU/MU' 'HP' 'LP' 'dBSPL' 'CenterRate' 'BehavState'
%'depth' 'jitter' 'baselineFR' 'DiscrmVar' 'trialN'
InfoTable.Session      = 'placeholder';
InfoTable.cluname      = 'placeholder';
InfoTable.UnitType     = discretize(4,[0.5 1.5 2.5 3.5 4.5],'categorical',{'unk', 'SU', 'MU', 'noise'});
InfoTable.HP           = 0;
InfoTable.LP           = 0;
InfoTable.dBSPL        = 0;
InfoTable.CenterRate   = 0;
InfoTable.BehavState   = categorical({'P'},{'P' 'D' 'A'},{'Passive' 'Drinking' 'Active Behavior'});
InfoTable.depth        = 0;
InfoTable.jitter       = 'placeholder';
InfoTable.baselineFR   = 0;
InfoTable.prev250msFR  = 0;
InfoTable.prev100msFR  = 0;
InfoTable.DiscrmVar    = 'placeholder';
InfoTable.trialN       = 0;


%% FOR EACH INDEPENDENT VARIABLE

for iV = { 'jitter' 'depth' }

    indVar = iV{1};
    nmfieldname = ['nm_' indVar];
    close all

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
                    
                    
                    %% Get spikes and save stim info
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    standardPd_collectdata( rasters,...
                        session, unit{:}, Data.(unit{:}).labels(1,2),...
                        stimpars, behavs(ib), indVar )
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    
                end %behavs
            end %pars
        end %stim block
    end %unit
end %session
end % for indVar

% Remove first row of Table if it's still the placeholder
if strcmp(InfoTable{1,1},'placeholder')
    InfoTable(1,:)=[];
    AmalgamatedSpikes(1,:)=[];
    OtherPdSpikes(1,:,:)=[];
end
AmalgamatedSpikes(size(InfoTable,1)+1:end,:)=[];
OtherPdSpikes(size(InfoTable,1)+1:end,:,:)=[];
if size(AmalgamatedSpikes,1) ~= size(InfoTable,1)
    keyboard
end

%% Save the data 
tablefn = sprintf('StandardPd_StimInfo_%s',subject);
writetable(InfoTable,fullfile(fn.processed,tablefn));

save(fullfile(fn.processed,'StandardPd_Spikes'),'AmalgamatedSpikes','-v7.3');
save(fullfile(fn.processed,'OtherPds_Spikes'),'OtherPdSpikes','-v7.3');


end





