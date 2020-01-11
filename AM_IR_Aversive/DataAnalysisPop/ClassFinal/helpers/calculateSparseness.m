function S = calculateSparseness(TuningCurve)

nStim = size(TuningCurve,1);

S = 1 - ( sum(TuningCurve./nStim)^2 ) / ( sum((TuningCurve.^2)./nStim) );

end

