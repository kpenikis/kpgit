function AssessUnits_VS(select_subject, select_session )
%
%  AssessUnits_VS( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   
%   Excludes datapoints based on: min Ntrials.
%
%   If no input arguments, runs through all units, overwriting previous
%   file. If a unit is specified in the input arguments, program adds the
%   unit to the existing files. 
%   If the unit entered as the input already exists, the program will exit. 
%   In the future, could allow the row to update individually.
%   
%  KP, 2018-06
%


SAVENAME = 'UnitsVS';


close all
global fn 

fn = set_paths_directories('','',1);

rng('shuffle')
%!!!!!!!!!!!!!!!!!
FRcutoff =  0.01;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!


%%  Prepare UnitInfo table and UnitData struct

UnitData = struct;

UnitInfo = table;
UnitInfo.Subject   = ' ';
UnitInfo.Session   = ' ';
UnitInfo.Channel   = nan;
UnitInfo.Clu       = nan;

N=0;

%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    if ~iscell(select_subject)
        subjects = {select_subject};
    end
else
    subjects = { 'AAB_265054'  };
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    

%%  SESSIONS

% Get list of sessions to check for sorted data

if nargin>1 && exist('select_session','var')
    if ~iscell(select_session)
        Sessions = {select_session};
    end
else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*VS_Spikes.mat'));
    
    Sessions = [];
    for ifn = 1:numel(SpkFns)
        Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
    end
    
end
% Sessions = flipud(Sessions);


% Step through each session
for sess = Sessions'
    
session = char(sess);

% Load data files
clear Spikes Clusters
fprintf('\n***************************************\nLoading %s sess %s...\n',subject,session)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));


% FOR NOW, ADD PLACEHOLDER FOR ARTIFACT TRIALS FOR SYNAPSE/KS DATA
if ~isfield(Info,'artifact')
        Info.artifact(64).trials = [];
end

 
%% Identify clusters to process, depending on whether sorting was done with UMS or KS

if exist('Spikes','var')                                     % >>> UMS <<<
    
    load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/SortingConfig/NN_A4x4_16/geometry_NN_A4x4_16.mat','kcoords','chanMap');
    
    
    Channels = find(Spikes.man_sort==1);
    for ich = 1:length(Channels)
        
        clus = Spikes.sorted(Channels(ich)).labels(Spikes.sorted(Channels(ich)).labels(:,2)==2,1)';
        if isempty(clus)
            continue
        else
            
            for iclu = 1:numel(clus)
                
                channel    = Channels(ich);
                clu        = clus(iclu);
                shank      = kcoords(chanMap==channel);
                spiketimes = unique(Spikes.sorted(Channels(ich)).spiketimes(Spikes.sorted(Channels(ich)).assigns==clu') * 1000);  %ms
                
                if ~isempty(spiketimes)
                    assess_this_unit
                end
                
                clear channel clu shank spiketimes
                
            end
            
        end
    end
    
    
elseif exist('Clusters','var')                                % >>> KS <<<
    
    if ~isfield(Clusters,'maxChTemp')
        Clusters = KS_to_Spikes(subject,session);
    end
    
    for iclu = 1:numel(Clusters)
        
        channel    = Clusters(iclu).maxChannel;
        clu        = Clusters(iclu).clusterID;
        shank      = Clusters(iclu).shank;
        spiketimes = unique(Clusters(iclu).spikeTimes * 1000)'; %ms
        
        if ~isempty(spiketimes)
            assess_this_unit_VS
        end
        
        clear channel clu shank spiketimes
        
    end
    
                                             % or else something's wrong...
elseif exist('Clusters','var') && exist('Spikes','var')
    keyboard
else
    keyboard
    Clusters = KS_to_Spikes(subject,session);
end % Spikes data type

    
end %session
end %subject

UnitInfo(1,:)=[];


%% Measure waveform shapes

[UnitInfo,UnitData] = assessWaveformShapes( SAVENAME, [SAVENAME '_WaveformShapes'], 1, UnitInfo, UnitData);


%% Save Unit files 

save(fullfile(fn.processed,SAVENAME),'UnitInfo','UnitData','-v7.3');



end %function











