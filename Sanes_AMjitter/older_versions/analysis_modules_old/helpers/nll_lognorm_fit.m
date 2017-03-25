function nll = nll_lognorm_fit(test_mu,test_sigma,test_stimuli,dp_observed)

if test_mu<=0 || test_sigma<=0
    nll=nan;
    sprintf('skipped mu=%i sigma=%i', test_mu, test_sigma)
    return
end

dp_computed = lognorm(test_mu,test_sigma,test_stimuli);

nll_thisStim=nan(1,size(test_stimuli,2));
for istim=1:numel(test_stimuli)
    nll_thisStim(istim) = -(dp_observed(istim) .* log(dp_computed(istim)) + (1-dp_observed(istim)).* log(1-dp_computed(istim)));
end

nll = sum(nll_thisStim);
end


