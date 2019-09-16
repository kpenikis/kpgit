function avgRate = getPrecRates(pdStart,twin,RateStream)
% timestamps don't exactly line up with TrialData onsets
% check and adjust in future
%
% average in log space
%

avgRate = 2^ mean( log2(RateStream( floor(pdStart) + (-twin:-1) )) );
% avgRate = mean( RateStream( floor(pdStart) + (-twin:-1) ) );

end