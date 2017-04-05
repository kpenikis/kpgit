function fn = set_paths_directories(subject,session)

% Add Processing and Analysis folders to path
addpath(genpath('DataAnalysis'), genpath('DataProcessing'))


% Set basic folders based on where runnuing from 

[~,computername] = system('hostname');

if strncmp(computername,'Regina',5)  
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % PC at 1012 rig
    fn.raw = 'G:\RawData';
    fn.processed  = 'G:\ProcessedData_temp';
   
else
    
    if exist('/Volumes/Seagate-1_KP','dir')
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % external harddrive connected to macbook
        fn.raw = '/Volumes/Seagate-1_KP/RawData';
        fn.processed  = '/Volumes/Seagate-1_KP/ProcessedData_temp';
        
    else
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % local macbook harddrive
        fn.raw       = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData';
        fn.processed = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        
    end
end

% Set subject and session directories
if nargin==2
    fn.sess_data   = fullfile(fn.processed,subject,session);
    fn.rasterplots = fullfile(fn.processed,subject,'^rasters',session);
    fn.anplots     = fullfile(fn.processed,subject,'^an_plots',session);
    fn.nmplots     = fullfile(fn.processed,subject,'^nm_plots',session);
end

end