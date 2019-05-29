function FF  = bootstrap_for_FF( TrialResps, nTrs )
% called by predObserved 
% KP 2019-05
% 

warning('off', 'MATLAB:rankDeficientMatrix')
rng('shuffle')
Iterations = 1000;

FFi = nan(1,Iterations);

for ii = 1:Iterations
    
    Resps = TrialResps(randi(length(TrialResps),nTrs));
    FFi(ii) = var(Resps)/mean(Resps);
    
end

FF = mean(FFi);

end
