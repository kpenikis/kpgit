function Clusters = KS_to_Spikes(SUBJECT,SESS_LABEL)
% 

% Load necessary data 
fn = set_paths_directories(SUBJECT,SESS_LABEL);

load(fullfile(fn.processed,SUBJECT,sprintf( '%s_sess-%s_Info',SUBJECT,SESS_LABEL)));

try
    load(fullfile(fn.sessdata,'sorting','KS_Output.mat'));
    clust_id = readNPY(fullfile(fn.sessdata,'sorting','spike_clusters.npy'));
    try
        clust_group = importdata(fullfile(fn.sessdata,'sorting','cluster_groups.csv'));
    catch
        disp('Session hasn''t yet been manually sorted.')
        keyboard
    end
catch
    keyboard
end

try
    load(rez.ops.chanMapFile,'chanMap','connected', 'xcoords', 'ycoords','kcoords');
catch % bc i moved the probe files
    pathparts = strsplit(rez.ops.chanMapFile,'/');
    if any(strcmp(pathparts,'ProcessedData'))
        pathparts(strcmp(pathparts,'ProcessedData')) = '';
    end
    rez.ops.chanMapFile = ['/' fullfile(pathparts{:})];
    load(rez.ops.chanMapFile,'chanMap','connected', 'xcoords', 'ycoords','kcoords');
    if ~exist('connected','var')
        keyboard
    end
end


% Create output struct
clusterLabels = cell(length(clust_group) - 1, 2);

for ii = 1:length(clust_group) - 1
    line = textscan(clust_group{ii+1},'%d %s');
    clusterLabels{ii, 1} = line{2};%label assigned to the cluster
    clusterLabels{ii, 2} = line{1};%cluster id 
end

goodClusterNames = strcmp([clusterLabels{:,1}],'good'); %only take good/SUs
goodSpikeIndices = ismember(int32(clust_id), [clusterLabels{goodClusterNames,2}]);
goodClusters = unique(clust_id(goodSpikeIndices));

connChs = find(connected);
Clusters = struct('spikeTimes',[],'clusterID',[],'maxChannel',[],'coordinates',[],'shank',[]);

for ii = 1:length(goodClusters)
    
    clusterOfInterest = goodClusters(ii);
    clusterSpikeIndices = ismember(clust_id, clusterOfInterest);
    clusterSpikeTimes = rez.st3(clusterSpikeIndices,1)./Info.fs; %seconds, not rounded
    
    %Find (channel) location of max waveform
    tmpSpike = clust_id==clusterOfInterest;
    origC = unique(rez.st3(tmpSpike,2));
    meanTemplate = mean(rez.Wraw(:,:,origC),3);
    [~,idx_connChs] = max(mean(abs(meanTemplate ),2));
    
    idx_allChs = connChs(idx_connChs);
    
    Clusters(ii).spikeTimes  = unique(clusterSpikeTimes);
    Clusters(ii).clusterID   = clusterOfInterest;
    Clusters(ii).maxChannel  = chanMap(idx_allChs);
    Clusters(ii).coordinates = [xcoords(idx_allChs) ycoords(idx_allChs)];
    Clusters(ii).shank       = kcoords(idx_allChs);
    Clusters(ii).maxChTemp   = mean(rez.Wraw(idx_connChs,:,origC),3)./1e6;
    Clusters(ii).Amplitudes  = rez.st3(tmpSpike,3);
end



% Save the output file
saveDataDir = fullfile(fn.processed,SUBJECT);
savename = sprintf('%s_sess-%s_Spikes',SUBJECT,SESS_LABEL);
save(fullfile( saveDataDir, savename),'Clusters','-v7.3');





end
