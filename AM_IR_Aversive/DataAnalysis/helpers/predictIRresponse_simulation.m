function [IR_Prediction_sim,IR_Pred_std_sim,pvals] = predictIRresponse_simulation(FRtrials)

global AMrates trMin
rng('shuffle')
Iterations = 1000;

if numel(FRtrials)==9
    FRtrials = FRtrials(1:8);
end
if size(FRtrials,1)>size(FRtrials,2)
    FRtrials = FRtrials';
end

% Simulate IR trials for within unit statistics
ntrs = cellfun(@length,FRtrials);
nt = min([ntrs(ntrs>0) 20]);
if nt<trMin, keyboard, end

Sim_mean = nan(Iterations,1);
Sim_std  = nan(Iterations,1);

p_tt  = nan(Iterations,numel(7:numel(FRtrials)));
p_rs  = nan(Iterations,numel(7:numel(FRtrials)));

for ii = 1:Iterations
    
    % ---- Simulate with Pdc ----
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
    
    % ---- Draw random Irr + Stats ----
    for iir = 7:numel(FRtrials)
        if isempty(FRtrials{iir}), continue, end
        
        [~,p_tt(ii,iir-6)] = ttest2(FR_sim,FRtrials{iir}(1:nt));
        p_rs(ii,iir-6)    = ranksum(FR_sim,FRtrials{iir}(1:nt));
        
%         q025 = prctile(FR_sim,2.5);
%         q975 = prctile(FR_sim,97.5);
%         if mean(FRtrials{iir}(1:nt))<q025 || mean(FRtrials{iir}(1:nt))>q975
%             
%         end
    end

end

IR_Prediction_sim = mean(Sim_mean);
IR_Pred_std_sim   = mean(Sim_std);

pvals = median(p_rs,1);

end


