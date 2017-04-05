function create_Data_struct(subject,session)


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
    
    
    % Go through each cluster of this channel and create a variable for
    % that cluster and a field in the Data struct with the same name
    for iu = 1:size(labels,1)
        
        varname = sprintf('ch%i_clu%i',ic,labels(iu,1));
        
        % Field labels is 3 element vector: [ clu#  SU/MU  channel# ] 
        eval(sprintf('%s.labels = [labels(iu,:) %i];',varname,ic))
        eval(sprintf('Data.%s.labels = [labels(iu,:) %i];',varname,ic))
        
        % Save cluster file
        filename = sprintf( '%s_sess-%s_%s',subject,session,varname);
        if ~exist(fn.sess_data,'dir')
            mkdir(fn.sess_data)
        end
        save(fullfile(fn.sess_data,filename),varname,'-v7.3');

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
    
    % Add raster to Data struct
    Data.(unit{:}).raster = raster;
    
end


% Save Data file
fprintf('\n##Saving Data struct to file.\n')
filename = sprintf( '%s_sess-%s_Data',subject,session); 
save(fullfile(fn.processed,subject,filename),'Data','-v7.3');




end


