function plotExUnitData(SUBJECT,SESS_LABEL,Clu1,Clu2)
% 

% Load necessary data 
fn = set_paths_directories(SUBJECT,SESS_LABEL);

load(fullfile( fn.processed,SUBJECT, sprintf( '%s_sess-%s_Info', SUBJECT,SESS_LABEL)));
load(fullfile( fn.processed,SUBJECT, sprintf('%s_sess-%s_Spikes',SUBJECT,SESS_LABEL)));

try
    load(fullfile(fn.sessdata,'sorting','KS_Output.mat'));
    clu_id   = readNPY(fullfile(fn.sessdata,'sorting','spike_clusters.npy'));
    spiketimes = readNPY(fullfile(fn.sessdata,'sorting','spike_times.npy')); % these are in samples, not seconds
catch
    keyboard
end


% Check root directory
fbinary_parts = strsplit(rez.ops.fbinary,'/');
if strcmp(fbinary_parts{4},'GDFS')
    rez.ops.root = fullfile(fn.processed,SUBJECT,SESS_LABEL,'sorting');
    rez.ops.fbinary = fullfile(fn.processed,SUBJECT,SESS_LABEL,'sorting',fbinary_parts{end});
    chanMapFile_parts = strsplit(rez.ops.chanMapFile,'/');
    rez.ops.chanMapFile = fullfile(fn.sorting,chanMapFile_parts{end-1},chanMapFile_parts{end});
end
load(rez.ops.chanMapFile,'chanMap','connected', 'xcoords', 'ycoords','kcoords');


gwfparams.dataDir = rez.ops.root;        % KiloSort/Phy output folder
gwfparams.fileName = fbinary_parts{end}; % .dat file containing the raw    [(1:end-4) '_filt.dat'] ; 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = rez.ops.NchanTOT;        % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 200;                     % Number of waveforms per unit to pull out
gwfparams.spikeTimes = spiketimes;       % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clu_id;        % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
gwfparams.theseUnits = [Clu1; Clu2];

wf = getWaveForms(gwfparams);
% wf.waveForms = [ Clu  x  iWaveform  x  Channel  x  sample ]


idx1 = [Clusters.clusterID]==Clu1;
idx2 = [Clusters.clusterID]==Clu2;
Clusters(idx1).maxChannel
Clusters(idx2).maxChannel

find(chanMap==Clusters(idx1).maxChannel)
find(chanMap==Clusters(idx2).maxChannel)

ChsToPlot = [37 38 39 40 41 42];

figure;
for ich = 1:numel(ChsToPlot)
    
    subplot(3,2,ich)
    theseWFs = permute(wf.waveForms(1,:,ich,:),[4 2 1 3]);
    plot(mean(theseWFs,2),'g')
    hold on
    theseWFs = permute(wf.waveForms(2,:,ich,:),[4 2 1 3]);
    plot(mean(theseWFs,2),'b')
    xlim([1 size(theseWFs,1)])
    set(gca,'Color','none')
    axis square
    
end



% 
% sp1   = round(Clusters([Clusters.clusterID]==Clu1).spikeTimes * 1000)';
% sp2   = round(Clusters([Clusters.clusterID]==Clu2).spikeTimes * 1000)';

% figure;
% subplot(1,2,1);
% histogram([-diff(sp1) diff(sp1)],[-xlimval:xlimval],'EdgeColor','none','FaceColor','c','FaceAlpha',1);
% set(gca,'xlim',[-xlimval xlimval],'Color','none')
% axis square
% title(['Clu #' num2str(Clu1)])
% 
% subplot(1,2,2);
% histogram([-diff(sp2) diff(sp2)],[-xlimval:xlimval],'EdgeColor','none','FaceColor','g','FaceAlpha',1);
% set(gca,'xlim',[-xlimval xlimval],'Color','none')
% axis square
% title(['Clu #' num2str(Clu2)])
% 





end
