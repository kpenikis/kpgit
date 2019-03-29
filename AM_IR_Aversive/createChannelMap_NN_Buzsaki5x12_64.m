% createChannelMap_NN_Buzsaki5x12_64()
% 
% Set channel mapping info for Buzsaki_5x12 probe. 
% 
% KP 2018-08

fn = set_paths_directories;

Shank1  = [61 60 63 62 58 59 56 57 53 52 55 54];
Shank2  = [49 48 51 50 18 17 20 22 25 26 24 21];
Shank3  = [32 29 30 28 31 27 5 1 6 4 3 2];
Shank4  = [11 10 8 7 12 14 15 16 44 45 46 47];
Shank5  = [40 41 42 43 39 38 37 36 33 64 34 35];
linpart = [13 19 9 23];

% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. 

chanMap = [Shank1 Shank2 Shank3 Shank4 Shank5 linpart];


% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(64, 1); 
rmv_chs = [64 linpart];
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

kcoords = [ones(size(Shank1)) 2*ones(size(Shank2)) 3*ones(size(Shank3)) 4*ones(size(Shank4)) 5*ones(size(Shank5)) 3*ones(size(linpart))];



% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

% single shank: 
x_shank = [  -1  -1  -1 -1 -1 -1    1    1   1     1    1    1 ]*10; %0 & 20 um
y_shank = [ 145 125 105 85 65 45   35   55  75    95  115  135 ];    % um

% all shanks:
xcoords = ([ repmat(x_shank,1,5)  [0 0 0 0] ]  + 200*kcoords - 190);
ycoords =  [ repmat(y_shank,1,5)  max(y_shank)+200*[1:size(linpart,2)] ];


% remake kcoords so bonus channels are sorted individually
kcoords = [ones(size(Shank1)) 2*ones(size(Shank2)) 3*ones(size(Shank3)) 4*ones(size(Shank4)) 5*ones(size(Shank5)) (5+(1:length(linpart)))];


% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map. 


save(fullfile('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/SortingConfig/Buzsaki5x12_64/geometry_Buzsaki5x12_64.mat'),...
    'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords')
 