function create_Data_struct(subject,session)
% create_Data_struct(subject,session)
%  First function to call after spike sorting. 
%  Creates a file called Data, which is in structure format. There is a
%  field for every SU or MU cluster in the session. 
%  Within the field for each cluster, is a reference of the cluster data,
%  a field called "labels", containing a 3 element vector: 
%       [ clu#  SU/MU  channel# ] 
%  Next, the program calls plot_rasters, and creates the raster plots if 
%  they don't already exist. 
%  The other field created is stimdata, and this contains all of the
%  various stimulus dimensions tested in this session, in a format
%  organized by block, then sound params (HP LP AMrate dB). For each
%  unique condition of these params, the raster entries are separated and
%  saved. The remaining params that could vary are jitter, depth, and
%  behavioral state. The rasters are saved in a field that reflects which
%  feature was varied and can be discriminated (jitter or depth, both only 
%  discriminated from 0 for now.)
%   
%  KP 2017-03
%

% Load Spikes file
fn = set_paths_directories(subject,session);
fprintf('loading Spikes file...\n')
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(fn.processed,subject,filename));


% Go through each channel and get SU and MU unit labels from Spikes file
Data=struct;
for ic = 1:numel(Spikes.sorted)
    
    % Get all single (2) and multi (3) unit cluster labels
    labels = [Spikes.sorted(ic).labels; 0 4];
    labels = labels((labels(:,2)==2)|(labels(:,2)==3),:);
    
    % Go through each cluster of this channel and create a field in the
    % Data struct with the same name
    for iu = 1:size(labels,1)
        % labels is 3 element vector: [ clu#  SU/MU  channel# ] 
        varname = sprintf('ch%i_clu%i',ic,labels(iu,1));
        eval(sprintf('Data.%s.labels = [labels(iu,:) %i];',varname,ic))
    end
    
end

% Add rasters to Data structure
allclusters = fieldnames(Data);
for unit = allclusters'
    
    channel = Data.(unit{:}).labels(1,3);
    clu     = Data.(unit{:}).labels(1,1);
    
    % Check if raster file already exists (just sessions LA and LC)
    clear raster
    filename = sprintf('%s_sess-%s_raster_ch%i_clu%i.mat',subject,session,channel,clu);
    foldername = sprintf('%s_sess-%s_rasters',subject,session);
    
    if (exist(fullfile(fn.processed,subject,foldername,filename),'file'))~=0
        load(fullfile(fn.processed,subject,foldername,filename))
        
    else
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         keyboard %pause to visually inspect plots
        close all
        raster = plot_rasters_wStim(subject,session,channel,clu);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
    
    
    %% Oraganize stimuli within raster
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    stimdata = organize_raster_stim(raster);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %% Add raster to Data struct
    Data.(unit{:}).stimdata = stimdata;
    
    
end



%% Save Data file
fprintf('\n##Saving Data struct to file.\n')
filename = sprintf( '%s_sess-%s_Data',subject,session); 
save(fullfile(fn.processed,subject,filename),'Data','-v7.3');




end


