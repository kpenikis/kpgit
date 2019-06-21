function [IR_Prediction_sim,IR_Pred_std_sim,pvals] = predictIRresponse_simulation(FRtrials)

global AMrates minTrs
Iterations = 200;

% Simulate IR trials for within unit statistics
ntrs = cellfun(@length,FRtrials);
nt = min(ntrs(ntrs>0));
if nt<minTrs, keyboard, end

Sim_mean = nan(Iterations,1);
Sim_std  = nan(Iterations,1);

for ii = 1:Iterations
    
    FR_sim = zeros(nt,5);
    for ir = 1:5
        these_trs = randperm(ntrs(ir+1),nt);
        for it = 1:nt
            tt = these_trs(it);
            FR_sim(it,ir) = FRtrials{ir+1}(tt) * ((1./AMrates(ir))/sum(1./AMrates));
        end
    end
    FR_sim = sum(FR_sim,2);
    
    Sim_mean(ii) = mean(FR_sim);
    Sim_std(ii)  = std(FR_sim);
end

IR_Prediction_sim = mean(Sim_mean);
IR_Pred_std_sim   = mean(Sim_std);


% ---- Stats ----
thisNIR = 0;
ptt   = nan(1,numel(7:numel(FRtrials)));
p_rs  = nan(1,numel(7:numel(FRtrials)));

for iir = 7:numel(FRtrials)
    if isempty(FRtrials{iir}), continue, end
    thisNIR = thisNIR+1;
    [~,ptt(thisNIR)] = ttest2(FR_sim,FRtrials{iir}(1:nt));
    p_rs(thisNIR)    = ranksum(FR_sim,FRtrials{iir}(1:nt));
end

pvals = p_rs;

end


