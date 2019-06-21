function fn = set_paths_directories(subject,session,UMS_flag)


% global SortingFlag

% Add Processing and Analysis folders to path
% cd /Users/kpenikis/Documents/kpgit/AM_IR_Aversive
addpath(genpath('DataAnalysis'),genpath('DataProcessing'))


% Set basic folders based on where runnuing from 

[~,computername] = system('hostname');

if strncmp(computername,'Regina',5)  
    keyboard
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % PC at 1012 rig
    basefolder  = 'G:\NYUDrive\Sanes\DATADIR\AM_IR_aversive'; %'G:\ProcessedData_temp';
    
    includePath='C:\gits\kpgit\ums2k_02_23_2012';
    
elseif strncmp(computername,'william',5)  
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % iMac in 1006
    basefolder  = '/Volumes/KPspace/GDFS/My Drive/Sanes/DATADIR/AMaversive';
    includePath='/Local/Users/sanesuser/Documents/MATLAB/KiloSort';
    
else
    if nargin==3 && UMS_flag
        includePath='/Users/kpenikis/Documents/kpgit/ums2k_02_23_2012';
%         removePath ='/Users/kpenikis/Documents/MATLAB/KiloSort';
    else
        includePath='/Users/kpenikis/Documents/MATLAB/KiloSort';
%         removePath ='/Users/kpenikis/Documents/kpgit/ums2k_02_23_2012';
    end
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    basefolder = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive';
%     basefolder = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/RawData/WWWlb_253396
    
    if exist(basefolder,'dir')==0
        keyboard
%         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         % local macbook harddrive
%         basefolder = '/Users/kpenikis/Documents/SanesLab/LocalData/AM_IR_aversive';        
    end
    addpath(genpath('/Users/kpenikis/Documents/MATLAB/Colormaps'))
end

% First, check if is in the current path (to avoid warning message)
% rmpath(genpath(removePath));

path(genpath(includePath),path);


% Set basic folders
fn.root       = basefolder;
fn.raw        = fullfile(basefolder,'RawData');
fn.processed  = fullfile(basefolder,'ProcessedData');
fn.stim       = fullfile(basefolder,'Stimuli');
fn.sorting    = fullfile(basefolder,'SortingConfig');
fn.figs       = fullfile(basefolder,'Figures');

% Set subject and session directories
if nargin==2
    fn.sessdata    = fullfile(fn.processed,subject,session);
    fn.rasterplots = fullfile(fn.processed,subject,'^rasters',session);
end

end