% createChannelMap_NN_A4x16_linear703()
% 
% Set channel mapping info for NN A4x16 linear probe. 
% 
% KP 2019-01

fn = set_paths_directories;

% bottom to top
Shank1  = [62 49 63 48 60 51 61 50 58 53 59 52 56 55 57 54];
Shank2  = [18 19 17 23 20 27 22 31 21 32 24 29 26 30 25 28];
Shank3  = [13 16  9 15  5 14  1 12  2 11  3 10  4  8  6  7];
Shank4  = [47 33 46 64 45 34 44 35 43 36 42 37 41 38 40 39];


% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel. 

chanMap = [Shank1 Shank2 Shank3 Shank4];


% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = true(64, 1); 
rmv_chs = 64;
for ii = 1:numel(rmv_chs)
    connected(chanMap==rmv_chs(ii)) = 0;
end


% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.

kcoords = [ones(size(Shank1)) 2*ones(size(Shank2)) 3*ones(size(Shank3)) 4*ones(size(Shank4))];



% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

% single shank: 
x_shank = zeros(size(Shank1)); 
y_shank = 50:50:800;    % um

% all shanks:
xcoords =  repmat(x_shank,1,4) + 500*(kcoords-1);
ycoords =  repmat(y_shank,1,4) ;



% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map. 


save(fullfile('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/SortingConfig/A4x16_linear703/geometry_A4x16_linear703.mat'),...
    'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords')
 