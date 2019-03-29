function CV = calc_binned_bootstrapped_CV(full_raster,binsize)


% Clip raster at closest integer multiple of binsize
raster = full_raster(:, 1:floor(size(full_raster,2)/binsize)*binsize );

% Bin raster, count spikes per bin for each trial
raster_bin = nan(size(raster,1),floor(size(full_raster,2)/binsize));
for it=1:size(raster,1)
    raster_bin(it,:) = sum(reshape(raster(it,:),binsize,size(raster,2)/binsize),1);
end

% Bootstrap 10 trials at a time, for a fair comparison of variance across stimuli
nTrs = min(10,size(raster,1));
nIterations = min(20*size(raster,1),500);
CVit = nan(nIterations,1);
for iteration = 1:nIterations
    
    trials = randperm(size(raster,1));
    
    CVit(iteration) = mean( std(raster_bin(trials(1:nTrs),:),1) ./ mean(raster_bin(trials(1:nTrs),:),1) ,'omitnan');
    
end


CV = mean(CVit);


end