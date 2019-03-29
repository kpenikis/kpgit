function [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,alpha)

% Bonferroni corrected significance

m = numel(pvals);

bonferroni_corrected_p = alpha/m;

iIR_sig_bonferroni = zeros(size(pvals));
iIR_sig_bonferroni(pvals<bonferroni_corrected_p) = 1;


% Bonferroni-Holm corrected significance

[pvals_sort,isort] = sort(pvals);

iIR_sig_bonferroniholm = zeros(size(pvals));
for k = 1:numel(pvals_sort)
    
    if pvals_sort(k) < (alpha / (m+1-k))
        iIR_sig_bonferroniholm(isort(k)) = 1;
    end
end


end

