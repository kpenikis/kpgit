function fn = set_paths_directories(subject,session)

% Add Processing and Analysis folders to path
addpath(genpath('DataAnalysis'), genpath('DataProcessing'))

% Set directory locations for loading and saving data
fn.raw       = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData';
fn.processed = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';

if nargin==2
    fn.sess_data   = fullfile(fn.processed,subject,session);
    fn.rasterplots = fullfile(fn.processed,subject,'^rasters',session);
    fn.anplots     = fullfile(fn.processed,subject,'^an_plots',session);
    fn.nmplots     = fullfile(fn.processed,subject,'^nm_plots',session);
end

end