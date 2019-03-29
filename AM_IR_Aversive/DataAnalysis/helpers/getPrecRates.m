function avgRate = getPrecRates(pdStart,allPds,ipd,pstid,twin)

global AMrates

avgRate = nan;
return 

if pdStart>0
    prevStim = [0 1000./allPds(1:ipd-1)];
    prevStim = cumsum([prevStim repmat(1000/AMrates(pstid-1),1,AMrates(pstid-1))]);
else
    if pstid<7
        prevStim = [0 cumsum(repmat(1000/AMrates(pstid-1),1,AMrates(pstid-1)))];
    else
        keyboard
    end
end
ii = find(prevStim>=twin,1);

avgRate = sum( 1000./diff(prevStim(1:ii)) .* [diff(prevStim(1:ii-1)) twin-sum(diff(prevStim(1:ii-1)))] ) / twin;

end