
function create_Data_struct(subject,session)


% Load Spikes file
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
fprintf('loading Spikes file...\n')
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(datadir,subject,filename));


% Go through each channel and get SU and MU unit labels from Spikes file
Data.ch=struct;
for ic = 1:numel(Spikes.sorted)
    
    % Get all single (2) and multi (3) unit cluster labels
    labels = [Spikes.sorted(ic).labels; 0 4];
    labels = labels((labels(:,2)==2)|(labels(:,2)==3),:);
    
    % Go through each cluster of this channel
    for iu = 1:size(labels,1)
        ch(ic).clu(iu).label = labels(iu,:);
    end
    
end

% Add to Data struct
Data.ch = ch;


% Now add raster struct in proper place
for channel = 1:numel(Data.ch)
    for iu = 1:numel(Data.ch(channel).clu)
        
        clu = Data.ch(channel).clu(iu).label(1);
                
        % Check for raster file
        
        filename = sprintf('%s_sess-%s_raster_ch%i_clu%i.mat',subject,session,channel,clu);

        if (exist(fullfile(datadir,subject,filename),'file'))~=0
            
            load(fullfile(datadir,subject,filename))
            
        else
            keyboard
            close all
            raster = pp_plot_rasters_wStim(subject,session,channel,clu);
            
        end
        
        
        % Add raseter to Data struct
        Data.ch(channel).clu(iu).raster = raster;
        clear raster
        
    end
    
end

% Save Data file
fprintf('\n##Saving Data struct to file.\n')
filename = sprintf( '%s_sess-%s_Data',subject,session); 
save(fullfile(datadir,subject,filename),'Data','-v7.3');



end


