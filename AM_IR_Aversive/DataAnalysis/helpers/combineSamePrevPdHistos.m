function FinalData = combineSamePrevPdHistos(MPHhistos,MPH)

global AMrates

% Combine rows that have the same previous period
unqPrPds = unique(MPH.PrevPd)';
FinalData = nan(5,size(MPHhistos,2));
for ipp = unqPrPds
    theseidx = MPH.PrevPd==ipp;
    FinalData(ipp==AMrates,:) = mean(MPHhistos(theseidx,:),1);
end

end