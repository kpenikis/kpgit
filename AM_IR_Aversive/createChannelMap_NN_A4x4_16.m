% createChannelMap_NN_A4x4_16()
% 
% Set channel mapping info for 16 ch probe. 
% 
% KP 2019-04

fn = set_paths_directories;

Shank1  = [9 10 4 12];
Shank2  = [3 2 11 1];
Shank3  = [7 6 8 14];
Shank4  = [15 16 13 5];

% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. 

chanMap = [Shank1 Shank2 Shank3 Shank4];


% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(16, 1); 
rmv_chs = 8;
for ii = 1:numel(rmv_chs)
    connected(chanMap==rmv_chs(ii)) = 0;
end


% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 

kcoords = [ones(size(Shank1)) 2*ones(size(Shank2)) 3*ones(size(Shank3)) 4*ones(size(Shank4))];



% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

% single shank: 
x_shank = [   0   0   0  0 ];
y_shank = [ 600 400 200  0 ];  % um

% all shanks:
xcoords =  repmat(x_shank,1,4)  + 200*kcoords - 200;
ycoords =  repmat(y_shank,1,4);



% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map. 


save(fullfile('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/SortingConfig/NN_A4x4_16/geometry_NN_A4x4_16.mat'),...
    'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords')
 