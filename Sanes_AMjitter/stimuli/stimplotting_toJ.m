

% raster is a structure with a field for each unique stimulus, containing
% information about the stim parameters. 
% with crucial variables: AMonset (ms) and AM depth

% call this withing plotting code...
Wave = ap_stimplotting(subject,raster);
Wcolor = [0.85 0.85 0.85];

% and this adds the stim waveform to the plot
fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*ymaxval + (ymaxval)/2 ,...
    Wcolor,'EdgeColor','none')

