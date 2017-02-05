function nll = nll_weibull_fit(test_alpha,test_beta,test_stimuli,dp_observed)

if test_alpha<=0 || test_beta<=0
    nll=nan;
    sprintf('skipped a=%i B=%i', test_alpha, test_beta)
    return
end

dp_computed = weibull(test_alpha,test_beta,test_stimuli,1);

nll_thisStim=nan(1,size(test_stimuli,2));
for istim=1:numel(test_stimuli)
    nll_thisStim(istim) = -(dp_observed(istim) .* log(dp_computed(istim)) + (1-dp_observed(istim)).* log(1-dp_computed(istim)));
end

nll = sum(nll_thisStim);
end