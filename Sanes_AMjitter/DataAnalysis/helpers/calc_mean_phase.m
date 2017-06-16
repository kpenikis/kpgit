function mean_phase = calc_mean_phase(spikes)

phases = find(sum_spikes) .* (360/size(spikes,2));

sum_spikes = sum(spikes,1);

weighting = sum_spikes(sum_spikes>0);

keyboard

mean(phases.*weighting)



end