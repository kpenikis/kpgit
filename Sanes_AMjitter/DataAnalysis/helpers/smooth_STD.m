function SDsmooth = smooth_STD(spikes,binsize)
% currently outputting standard error of the mean!!

SDsmooth = nan(1,size(spikes,2)-binsize+1);

for ims = 1:(size(spikes,2)-binsize+1)
    
    t = ims:ims+binsize-1;
    SDsmooth(ims) = std(mean(spikes(:,t),2)*1000);
    
end

% Fill out to full duration (when spikes has no buffer)
buffer = (size(spikes,2)-length(SDsmooth))/2;
SDsmooth = [repmat(SDsmooth(1),1,buffer) SDsmooth repmat(SDsmooth(end),1,buffer)]/sqrt(size(spikes,1));


end



