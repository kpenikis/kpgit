function make_Atlas_ChannelMap(fpath)
% create a channel Map file for simulated data (eMouse)

% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
% be dead channels. 

chanMap = 1:16;

% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(16, 1); connected(16) = 0;

% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

xcoords = 2 * 1+[ 16 14 15 10 11 10 9 6 5 4 5 0 1 -1 0 nan];
ycoords = 1 *   [  2  2  1  3  2  1 2 2 3 2 1 3 2  2 1 nan]; 

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords = [ 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 nan];

% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map. 

save(fullfile(fpath, 'geometry_Atlas4tetprobe.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords')