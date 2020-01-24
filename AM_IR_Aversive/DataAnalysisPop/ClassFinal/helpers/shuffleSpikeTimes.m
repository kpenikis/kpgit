function CTTS = shuffleSpikeTimes(CTTS)
% shuffles order of time dimension

for i1 = 1:size(CTTS,1)
    for i3 = 1:size(CTTS,3)
        for i4 = 1:size(CTTS,4)
            for i5 = 1:size(CTTS,5)
                CTTS(i1,:,i3,i4,i5) = CTTS(i1,randperm(size(CTTS,2)),i3,i4,i5);
            end
        end
    end
end



