function [pkFRsort,ipkFR] = rankPeakFR(CTTS)


% Get peak FR for each stimulus 
stim_pks=nan(size(CTTS,1),8);
for iUn = 1:size(CTTS,1)
    stim_pks(iUn,:) = permute( 1000.*max(mean(CTTS(iUn,:,:,:),3,'omitnan'),[],2) ,[3 4 1 2]);
end

% Rank by median of peak FRs across all stimuli 
[pkFRsort,ipkFR] = sort(median(stim_pks,2),'descend');

if any(isnan(pkFRsort))
    keyboard
end
% Rank by mean of peak FRs across all stimuli 
% [pkFRsort,ipkFR] = sort(mean(stim_pks,2),'descend');


end