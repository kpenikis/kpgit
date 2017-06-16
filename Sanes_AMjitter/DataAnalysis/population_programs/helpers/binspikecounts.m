function spks_binned = binspikecounts(data,binsize)
% Data input must be trials X time in ms.

spks_binned = nan(size(data,1),size(data,2)/binsize);

% Step through each row, corresponding to a trial or a unit
for it = 1:size(data,1)
    spks_binned(it,:) = sum(reshape(data(it,:),[binsize size(data(it,:),2)/binsize]),1);
end

end