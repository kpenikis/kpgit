
% compare_std_range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Things to note: 
%  - distributions are based on period, not rate, so are assymmetric
%      (the resulting mean rate is NOT the desired rate, but a smidge higher!)
%  - the distributions are clipped according to fraction of the period, and
%    NOT percentile of distribution, which is why the range asymptotes.
%    By cutting at 0.25*(1/AM_Hz) is the 7th and 93rd percentile of the
%    distribution for std=1/6*(1/AM_Hz).
%  - checked and the amount that std is reduced IS consistent for each rate
%    used in the Overath paper (3 29 57 Hz). effective std is closer to
%    12.5% than 16.7%.
%  


% make and plot ex rate distributions with tobias's methods
AMrates = [4]; %10 15 20 30 60 100 

for ir = 1:numel(AMrates)
    
AM_Hz = AMrates(ir);
dur = 6; %s

% set_stds = [1/4 0.2 1/6 1/8 0.1 1/12 1/16 1/20 1/25 1/50 0];
% set_stds = [0:0.01:0.25]; %
set_stds = 1/6;

clear data
data = struct;

for isd = 1:numel(set_stds)
std_Hz = set_stds(isd);

reps=100;
ranges = nan(reps,2);
for irep = 1:reps

AM_jitt = randn(1,100*dur*AM_Hz)*(1/AM_Hz)*std_Hz +1/AM_Hz;  % create normal distribution with std = 1/6 and mean = 1/AM_Hz; 100 times more values to choose from to help with selecting
AM_jitt = AM_jitt - mean(AM_jitt) + 1/AM_Hz;    % make sure distribution has mean of 1/AM_Hz because of sampling error
  AM_jitts = AM_jitt(find(AM_jitt < 1/AM_Hz+.25*(1/AM_Hz)));  % only include AM rates below 1/AM_Hz+.25*(1/AM_Hz)
  AM_jitts = AM_jitts(find(AM_jitts > 1/AM_Hz-.25*(1/AM_Hz)));    % only include AM rates above 1/AM_Hz-.25*(1/AM_Hz)
% AM_jitts = AM_jitt((1./AM_jitt)>0);
AM_select = randperm(length(AM_jitts)); % additional randomization step (not really necessary)
AM_jitter = AM_jitts(AM_select(1:round(dur*AM_Hz)));


[counts,pds] = hist(AM_jitt,100);
[counts_clipped,pds_clipped] = hist(AM_jitts,100);

% figure(1); clf; hold on
% plot(1./pds, counts./max(counts),'b','LineWidth',1)
% plot(1./pds_clipped, counts_clipped./max(counts_clipped),'k','LineWidth',1)
% xlim([AM_Hz-2 AM_Hz+8])

% fprintf('\nAM rate %i Hz, fraction full: %2.5f, fraction clipped: %2.5f',AM_Hz,std(AM_jitt)/(1/AM_Hz),std(AM_jitts)/(1/AM_Hz))

% Calculate range width 
ranges(irep,:) = 1./[max(pds_clipped) min(pds_clipped)];

% pause(1)
end
end
end



% figure;
% yranges = [1:reps;1:reps];
% plot(ranges',yranges)
% xlim([AM_Hz-1 AM_Hz+2])
% ylim([0 51])
% set(gca,'XScale','log')


data(isd).set_std = std_Hz;
data(isd).range   = mean(diff(ranges'));
data(isd).rng_std = std(diff(ranges'));


end


figure(6); hold on
plot([data.set_std],[data.range],'o-r')
xlim([0 max(set_stds)])
ylim([-0.1 max([data.range])*1.2])
xlabel('std (prop. of AM pd) from tobias stimulus set')
ylabel('corresponding range of AM rates (Hz)')
title(['std - range relationship for ' num2str(AM_Hz) ' Hz AM rate'])



end



