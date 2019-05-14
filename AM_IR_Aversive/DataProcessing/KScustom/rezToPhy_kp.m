function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy_kp(rez)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%

% rez.ops.chanMapName = ops.chanMap;
% [length(unique(rez.st3(:,2))),length(rez.st3(:,2))]
savePath=fullfile(rez.ops.root);
% spikeTimes will be in samples, not seconds

% if 0%check for doubled spike times, this is ugly and old
%     for i=unique(rez.st3(:,2))'
%         inds=rez.st3(:,2)==i;
%         doubled(i)=length(rez.st3(inds,1))-length(unique(rez.st3(inds,1)));
%     end
%     if sum(doubled>10)
%         doubled(logical(doubled))
%     end
% end

%If you already ran phy on your data, it will have a bad config file.
%delete it. 
fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
   delete(fullfile(savePath, fs(i).name)); 
end
fs = dir(fullfile(savePath, '*.csv'));
if ~isempty(fs)
    delete(fullfile(savePath, fs.name)); 
end
fs = dir(fullfile(savePath, '*.log'));
if ~isempty(fs)
    delete(fullfile(savePath, fs.name)); 
end
dat_folder=fileparts(savePath);
if exist(fullfile(dat_folder,'.phy'),'dir')
    rmdir(fullfile(dat_folder,'.phy'),'s')
end


spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4 %what is this for?
    spikeClusters = uint32(1+rez.st3(:,5));
end
clusterIDs=[];
amplitudes = rez.st3(:,3);
mV=getmV(rez);
Nchan = rez.ops.Nchan;

% rez.connected(8) = 2;
connected   = rez.connected(:);
xcoords     = rez.xcoords(:); %for display in template_gui
ycoords     = rez.ycoords(:); %for display in template_gui
kcoords     = rez.ops.kcoords(:); 
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
for iNN = 1:rez.ops.Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

%%
if ~isempty(savePath)
    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(mV, fullfile(savePath, 'amplitudes_mV.npy'));

    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
%     Fs = rez.ops.fs;
    conn        = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savePath, 'channel_map.npy'));
%     writeNPY(int16(kcoords(conn)), fullfile(savePath, 'channel_shank_map.npy'));
    writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savePath, 'channel_positions.npy'));
%     writeNPY([xcoords ycoords], fullfile(savePath, 'channel_positions.npy'));

    writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    whiteningMatrix = rez.Wrot/200;
    whiteningMatrixInv = whiteningMatrix^-1;
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end
    
    
    % If a filtered file was created, save path and point phy there instead
%     if isfield(rez.ops,'fbinary_filt')
%         datpath = rez.ops.fbinary_filt;
%     else 
        datpath = rez.ops.fbinary;
%     end
    
    
    % Make params file
    fid = fopen(fullfile(savePath,'params.py'), 'w');
    fprintf(fid,['dat_path = ''' datpath '''\n']);
    fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
    fprintf(fid,'dtype = ''int16''\n');
    fprintf(fid,'offset = 0\n');
    if mod(rez.ops.fs,1)
        fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
    else
        fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
    end
    fprintf(fid,'hp_filtered = False');
    fclose(fid);
end
disp('done writing to phy')


%%
figure(1);clf;
tS=rez.st3(:,1)/rez.ops.fs/60;
[n,edges]=histcounts(tS,linspace(0,max(tS),100));
dT=diff(edges(1:2));
n=n/length(unique(rez.st3(:,2)))/60/dT;
plot(edges(2:end-2),n(2:end-1));
axis tight;
xlabel('minutes')
ylabel('mean FR (Hz)')


%%
s=whos('rez');
if (s.bytes/1e9)<2
    save(fullfile(rez.ops.root,'KS_Output.mat'),'rez');
else
    save(fullfile(rez.ops.root,'KS_Output.mat'),'rez','-v7.3');
end
