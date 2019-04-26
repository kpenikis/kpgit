function avgRate = getPrecRates(pdStart,twin,RateStream)
% timestamps don't exactly line up with TrialData onsets
% check and adjust in future
% 

avgRate = mean( RateStream( floor(pdStart) + (-twin:-1) ) );

end