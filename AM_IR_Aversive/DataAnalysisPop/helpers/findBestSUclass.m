
% Find best SU classifier


load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/StimClass/Pdc/PCMat_15trTemp.mat')

PCs = nan(size(PCMat,3),size(PCMat,1));

for iUn = 1:size(PCMat,3)
    
    PCs(iUn,:) = diag(PCMat(:,:,iUn));
    
end


[maxMean,iMaxMean]     = max(mean(PCs,2,'omitnan'))
[maxMedian,iMaxMedian] = max(median(PCs,2,'omitnan'))


% unit 40
% unit 115



