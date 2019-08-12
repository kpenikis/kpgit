
% code to get raw waveforms


mmf = memmapfile(params.filename, 'Format', {params.dataType, [nChInFile nSamp], 'x'});

st = readNPY('spike_times.npy'); % these are in samples, not seconds
clu = readNPY('spike_clusters.npy');

theseST = st(clu==19); % spike times for cluster 19

extractST = theseST(1:min(100,length(theseST))); %extract at most the first 100 spikes

nWFsToLoad = length(extractST);

nCh = 64; % number of channels
wfWin = [-30:30]; % samples around the spike times to load
nWFsamps = length(wfWin);

theseWF = zeros(nWFsToLoad, nCh, nWFsamps);

for i=1:nWFsToLoad
    tempWF = mmf.Data.x(1:nChInFile,extractST(i)+wfWin(1):extractST(i)+wfWin(end));
    theseWF(i,:,:) = tempWF(params.chanMap+1,:);
end

