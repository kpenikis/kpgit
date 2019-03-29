function AssessUnits(select_subject, select_session, select_channel, select_clu)
%
%  AssessUnits( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%
%   If no input arguments, runs through all units, overwriting previous
%   file. If a unit is specified in the input arguments, program adds the
%   unit to the existing files. 
%   If the unit entered as the input already exists, the program will exit. 
%   In the future, could allow the row to update individually.
%   
%  KP, 2018-06
%

close all
global fn 


SAVENAME = 'Units_AAB055';


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)
scrsz = get(0,'ScreenSize');

hf=figure; hold on
hs_spk = subplot(1,3,1);
ylim([0 4*360])
axis square
hs_gap = subplot(1,3,2);
ylim([0 4*360])
axis square
hs_diff= subplot(1,3,3); hold on
axis square
plot([0 34],[0 0],'--k')
hold off

set(hf,'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])


fn = set_paths_directories('','',1);

%  Load SplitUnits
load([fn.processed '/SplitUnits.mat'])


rng('shuffle')
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  5;
%!!!!!!!!!!!!!!!!!


%%  Prepare UnitInfo table and UnitData struct

if nargin<1 %start over, re-running all units
    
    RERUN = input('If want to re-run ALL UNITS, enter 1. Otherwise, enter 0.   ');
    
end

if RERUN == 1
    
    UnitData = struct;
    
    UnitInfo = table;
    UnitInfo.Subject   = ' ';
    UnitInfo.Session   = ' ';
    UnitInfo.Channel   = nan;
    UnitInfo.Clu       = nan;
    UnitInfo.RespType  = ' ';    
    
    N=0;
    
else  %to add units to pre-existing struct
    
    if SUonly==1
        filename = 'Units';
    elseif SUonly==0
        filename = 'UnitInfo_MUtoo';
    end
    
    q = load(fullfile(fn.processed,filename));
    
    UnitData = q.UnitData;
    UnitInfo = q.UnitInfo;
    
    N = numel(UnitData);
    
end

allUn_FR_raw = nan(0,8);
allUn_FR_nrm = nan(0,8);
allUn_GR_raw = nan(0,8);


%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    if ~iscell(select_subject)
        subjects = {select_subject};
    end
else
    subjects = { 'AAB_265055' };%'AAB_265059' 'WWWf_253400' 'WWWlf_253395'  };
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    
    switch subject
        case 'WWWf_253400'
            subjcol = 'k';
        case 'WWWlf_253395'
            subjcol = 'k';
        case 'AAB_265059'
            subjcol = 'k';
    end
    

%%  SESSIONS

% Get list of sessions to check for sorted data

if nargin>1 && exist('select_session','var')
    if ~iscell(select_session)
        Sessions = {select_session};
    end
else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
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
fprintf('\n******************************\nLoading sess %s...\n',session)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));


% FOR NOW, ADD PLACEHOLDER FOR ARTIFACT TRIALS FOR SYNAPSE/KS DATA
if ~isfield(Info,'artifact')
        Info.artifact(64).trials = [];
end

 
%% Identify clusters to process, depending on whether sorting was done with UMS or KS

if exist('Spikes','var') && ~exist('Clusters','var')         % >>> UMS <<<
    
    Channels = find(Spikes.man_sort==1);
    for ich = 1:length(Channels)
        
        clus = Spikes.sorted(Channels(ich)).labels(Spikes.sorted(Channels(ich)).labels(:,2)==2,1)';
        if isempty(clus)
            continue
        else
            
            for iclu = 1:numel(clus)
                
                channel = Channels(ich);
                clu = clus(iclu);
                spiketimes = round(Spikes.sorted(Channels(ich)).spiketimes(Spikes.sorted(Channels(ich)).assigns==clu') * 1000);  %ms
                
                if ~isempty(spiketimes)
                    assess_this_unit
                end
                
                clear channel clu spiketimes
                
            end
            
        end
    end
    
    
elseif exist('Clusters','var') && ~exist('Spikes','var')      % >>> KS <<<
    
    Clusters = KS_to_Spikes(subject,session);
    
    for iclu = 1:numel(Clusters)
        
        channel    = Clusters(iclu).maxChannel;
        clu        = Clusters(iclu).clusterID;
        spiketimes = round(Clusters(iclu).spikeTimes * 1000)'; %ms
        
        if ~isempty(spiketimes)
            assess_this_unit
        end
        
        clear channel clu spiketimes
        
    end
    
    
                                             % or else something's wrong...
elseif exist('Clusters','var') && exist('Spikes','var')
    keyboard
else
    keyboard
    Clusters = KS_to_Spikes(subject,session);
    
end

    
end %session
end %subject



%% Finish and save UnitInfo and UnitData

if SUonly==1
    filename = 'Units';
elseif SUonly==0
    filename = 'Units_MUtoo';
end

if RERUN == 1
    UnitInfo(1,:)=[];
end


%% Combine data from split units and save files

combineSplitUnits( UnitData, UnitInfo, SAVENAME )  % ALSO SAVES UNITS FILE



return




% Label each unit as sparse or sustained (old)

separateSparseSustainedUnits(UnitData)


return






%% Play with data if wanted 

allUnResp_NRMraw = allUn_FR_raw./repmat(max(allUn_FR_raw,[],2),1,size(allUn_FR_raw,2));
figure;
imagesc(allUnResp_NRMraw)
colormap('bone')

figure;
imagesc(allUn_FR_nrm)
colormap('bone')


% Sort according to some feature

sortBy = 'pdc_raw_resp_range';
switch sortBy
    case 'warn_raw'
        [warn_raw,idx]=sort(Resp.raw(:,1),1);
    case 'raw_32'
        [raw_32,idx]=sort(Resp.raw(:,6),1);
    case 'nrm_32'
        [nrm_32,idx]=sort(Resp.nrm(:,6),1);
    case 'IR_pred'
        [IR_pred,idx]=sort(Resp.IR_pred,1);
    case 'diff_nrm_2_32'
        [diff_nrm_2_32,idx]=sort(Resp.nrm(:,2)-Resp.nrm(:,6),1);
    case 'pdc_raw_resp_range'
        [pdc_raw_resp_range,idx]=sort(range(Resp.raw(:,2:6),2),1);
    case 'pdc_nrm_resp_range'
        [pdc_nrm_resp_range,idx]=sort(range(Resp.nrm(:,2:6),2),1);        
end


%% Plots
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4) scrsz(3)/8 scrsz(4)];
figsize2 = [scrsz(3)/8 scrsz(4) scrsz(3)/8 scrsz(4)];

hf1 = figure;
set(hf1,'Position',figsize1,'NextPlot','add')
imagesc(allUn_FR_raw(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('mean FR resp (Hz)\nsorted by: %s',sortBy))

hf2 = figure;
set(hf2,'Position',figsize2,'NextPlot','add')
imagesc(allUn_FR_nrm(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('norm. mean resp, (z-FR)\nsorted by: %s',sortBy))


%% Identify units to exclude

[Resp(idx,1:4) table(eval(sortBy))]
figure; hist(abs(eval(sortBy)),10)



end %function










