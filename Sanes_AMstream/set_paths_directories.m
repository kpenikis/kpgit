function fn = set_paths_directories(subject,session)

% Add Processing and Analysis folders to path
addpath(genpath('DataAnalysis'), genpath('DataProcessing'))


% Set basic folders based on where runnuing from 

[~,computername] = system('hostname');

if strncmp(computername,'Regina',5)  
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % PC at 1012 rig
    basefolder  = 'G:\NYUDrive\Sanes\DATADIR\AMStream'; %'G:\ProcessedData_temp';
    
    includePath='C:\gits\kpgit\ums2k_02_23_2012';
    
else
    includePath='/Users/kpenikis/Documents/kpgit/ums2k_02_23_2012';
    
    if exist('/Volumes/KP_workingDataStore','dir')~=0
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % external harddrive connected to macbook
        basefolder = '/Volumes/KP_workingDataStore/NYUDrive/Sanes/DATADIR/AMStream'; %'/Volumes/Seagate-1_KP/RawData';
        
    else
        disp('External disk not detected; setting directories to macbook built in harddrive')
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % local macbook harddrive
        basefolder = '/Users/kpenikis/Documents/SanesLab/LocalData/AMStream';
        
    end
end
addpath(genpath(includePath));

% Set basic folders
fn.raw        = fullfile(basefolder,'RawData');
fn.processed  = fullfile(basefolder,'ProcessedData');
fn.stim       = fullfile(basefolder,'Stimuli');
fn.population = fullfile(fn.processed,'^Population');
fn.standardPd = fullfile(fn.population,'StandardPd');

% Set subject and session directories
if nargin==2
    fn.sess_data   = fullfile(fn.processed,subject,session);
    fn.rasterplots = fullfile(fn.processed,subject,'^rasters',session);
    fn.anplots     = fullfile(fn.processed,subject,'^an_plots',session);
    fn.nmplots     = fullfile(fn.processed,subject,'^nm_plots',session);
end

end