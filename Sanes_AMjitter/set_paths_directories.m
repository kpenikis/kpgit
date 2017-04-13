function fn = set_paths_directories(subject,session)

% Add Processing and Analysis folders to path
addpath(genpath('DataAnalysis'), genpath('DataProcessing'))


% Set basic folders based on where runnuing from 

[~,computername] = system('hostname');

if strncmp(computername,'Regina',5)  
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % PC at 1012 rig
    fn.raw = 'G:\NYUDrive\Sanes\DATADIR\AMJitter\RawData'; %'G:\RawData';
    fn.processed  = 'G:\NYUDrive\Sanes\DATADIR\AMJitter\ProcessedData'; %'G:\ProcessedData_temp';
    
    includePath='C:\gits\kpgit\ums2k_02_23_2012';
    
else
    includePath='/Users/kpenikis/Documents/kpgit/ums2k_02_23_2012';
    
    if exist('/Volumes/KP_workingDataStore','dir')~=0
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % external harddrive connected to macbook
        fn.raw = '/Volumes/KP_workingDataStore/NYUDrive/Sanes/DATADIR/AMJitter/RawData'; %'/Volumes/Seagate-1_KP/RawData';
        fn.processed  = '/Volumes/KP_workingDataStore/NYUDrive/Sanes/DATADIR/AMJitter/ProcessedData'; %'/Volumes/Seagate-1_KP/ProcessedData_temp';
        
    else
        disp('External disk not detected; setting directories to macbook built in harddrive')
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % local macbook harddrive
        fn.raw       = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData';
        fn.processed = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        
    end
end
addpath(genpath(includePath));

% Set subject and session directories
if nargin==2
    fn.sess_data   = fullfile(fn.processed,subject,session);
    fn.rasterplots = fullfile(fn.processed,subject,'^rasters',session);
    fn.anplots     = fullfile(fn.processed,subject,'^an_plots',session);
    fn.nmplots     = fullfile(fn.processed,subject,'^nm_plots',session);
end

end