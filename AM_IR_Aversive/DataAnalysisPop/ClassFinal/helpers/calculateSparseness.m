function S = calculateSparseness(TuningCurve)
% Input data should be [Nx1]


if size(TuningCurve,1)<size(TuningCurve,2)
    TuningCurve = TuningCurve';
end

nStim = size(TuningCurve,1);

S = 1 - ( sum(TuningCurve./nStim)^2 ) / ( sum((TuningCurve.^2)./nStim) );

end

