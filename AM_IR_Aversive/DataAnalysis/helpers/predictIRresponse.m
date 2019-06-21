function [Pred_FR,Pred_SEM] = predictIRresponse(PdcMeans,PdcVars,Ntrials)

global AMrates

Pred_FR     = sum( PdcMeans .* ((1./AMrates)/sum(1./AMrates)) );

% Errorbar option #1: weighted sum of SEMs
Pred_SEM    = sum( ( sqrt(PdcVars)./sqrt(Ntrials) ) .* ((1./AMrates)/sum(1./AMrates)) );

% Errorbar option #2: weighted sum of variances, then converted to overall SEM [--> var(a+b)=var(a)+var(b) ]
IR_Pred_sem2   = sqrt( sum( PdcVars .* ((1./AMrates)/sum(1./AMrates)) ) ) / sqrt(length(AMrates));



end