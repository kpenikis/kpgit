function Mdist = getPopMahalDistVector(Clusters,TrialData,spkshift,gausswinsize)
% Mdist = getPopMahalDistVector(Clusters,TrialData,gausswinsize)
% 
% First creates a matrix of activity in the population, smoothed with a 
% gaussian window. Then calculates the Mahalanobis distance of the vector 
% of activity in each ms, from the vectors of the rest of the session. 
% 
% Called by SessionRastersMahal.
% 
% KP, 2019-08
% 

FullSpikeMatrix = zeros(numel(Clusters),round(max(vertcat(Clusters.spikeTimes)*1000)));

for iClu = 1:numel(Clusters)
    
    % Get spiketimes (KS)
    spiketimes = unique(round(Clusters(iClu).spikeTimes*1000 - spkshift)');
    
    % Make gaussian smoothed version of activity
    Stream_FRsmooth = convertSpiketimesToFR(spiketimes,...
        length(FullSpikeMatrix),TrialData.onset(1),TrialData.offset(1),20,gausswinsize,'silence');
    
    % Add data to matrix
    FullSpikeMatrix(iClu,:)          = Stream_FRsmooth;
    
end

ms1 = TrialData.onset(1);
FullSpikeMatrix = FullSpikeMatrix(:,ms1:end);


Mdist = mahal(FullSpikeMatrix',FullSpikeMatrix');

Mdist = [nan(1,ms1-1) Mdist'];


end