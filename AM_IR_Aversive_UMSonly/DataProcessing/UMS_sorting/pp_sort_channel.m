function spikes = pp_sort_channel(data, thresh, reject)
%
%  spks = pp_sort_channel(data, thresh, reject)
%    Runs the UMS algorithm for data of a given channel/group of channels.
%    Important: include the optional inputs to ss_default_params, in order
%    to use the thresholds calculated from clean data segments.
% 
%  KP, 2016-04; last updated 2016-04
% 

global fs



% Run UMS2000 spike sorting algorithm

spikes = ss_default_params(fs,  'detect_method','manual',  'thresh',thresh,  'reject',reject);  %edit any default params?


% Make data struct into cell
for it = 1:size(data,1)
    data_cell{it} = data(it,:)';
end

spikes = ss_detect(data_cell,spikes);


% Run the rest of the sorting algorithm
try
    spikes = ss_align(spikes);
    spikes = ss_kmeans(spikes);
    spikes = ss_energy(spikes);
    spikes = ss_aggregate(spikes);
    
catch
    keyboard
end




end