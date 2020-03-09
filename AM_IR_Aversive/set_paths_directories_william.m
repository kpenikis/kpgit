function fn = set_paths_directories_william(subject,session,UMS_flag)

% Add Processing and Analysis folders to path
addpath(genpath('DataProcessing'),genpath('DataAnalysis'),genpath('DataAnalysisPop'))


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
    addpath(genpath('/Users/kpenikis/Documents/MATLAB/CircStat2012a'))
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



% Set some global variables that remain constant 

global trMin AMrates rateVec_AC rateVec_DB

trMin   = 10;
AMrates = [2 4 8 16 32];
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

end