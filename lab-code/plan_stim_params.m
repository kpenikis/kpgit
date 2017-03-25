function [ trialDur, experimentDuration_m, ml_total ] = plan_stim_params( n_jitter_trials, n_perdic_trials, ml_per_m )

if nargin<3
    n_jitter_trials = 20;
    n_perdic_trials = 10;
    ml_per_m = 0.27;
end

AMrates = 2.^[1:0.25:5];

nreps = 2;
trialDur = nreps*sum(1./AMrates);
totalDur_jitter = n_jitter_trials * trialDur; %sec

totalDur_perdic = ceil(length(AMrates)/4) *trialDur * n_perdic_trials;

ITI = 0;

experimentDuration_s = totalDur_jitter + totalDur_perdic + (ITI*n_jitter_trials + ITI*n_perdic_trials*ceil(length(AMrates)/4));
experimentDuration_m = experimentDuration_s/60;

ml_total = experimentDuration_m * ml_per_m;



figure; hold on
set(gca,'XScale','log')
plot(AMrates,2.*ones(size(AMrates)),'.b','MarkerSize',30)
plot(AMrates(1:2:end),3.*ones(size(AMrates(1:4:end))),'.k','MarkerSize',30)
ylim([0 10])
set(gca,'xtick',[1 2 4 8 16 32 64])



% rng('shuffle')
% rateVec = AMrates(randperm(length(AMrates)));

