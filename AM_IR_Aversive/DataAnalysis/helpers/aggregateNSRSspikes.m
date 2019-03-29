function [spiketimes_NS, spiketimes_RS] = aggregateNSRSspikes(UnitData,UnitInfo,Clusters)
% [spiketimes_NS, spiketimes_RS] = aggregateNSRSspikes(UnitData,UnitInfo,Clusters)
%
%

iNS = find(strcmp(UnitInfo.SpkShape,'NS'))';
iRS = find(strcmp(UnitInfo.SpkShape,'RS'))';

spiketimes_NS = [];
spiketimes_RS = [];

for iUn = 1:size(UnitInfo,1)
    
    % Get spiketimes (KS)
    iClu = find([Clusters.maxChannel] == UnitData(iUn).Channel(1) & [Clusters.clusterID] == UnitData(iUn).Clu(1));
    
    if ismember(iUn,iNS)
        spiketimes_NS = [spiketimes_NS (Clusters(iClu).spikeTimes*1000)'];
    elseif ismember(iUn,iRS)
        spiketimes_RS = [spiketimes_RS (Clusters(iClu).spikeTimes*1000)'];
    end
end

spiketimes_NS = sort(spiketimes_NS);
spiketimes_RS = sort(spiketimes_RS);

end