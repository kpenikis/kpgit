function Spikes = pp_sort_session( subject, session )
% 
%  pp_sort_session( subject, session) 
%    Loads Phys struct to run data through UMS spike sorting algorithm, and
%    creates the main output Spikes.clustered, which contains all relevant 
%    sorting info and results from UMS. Note: Spikes.sorted is created, by
%    replicating Spikes.clustered, and this will be where manual sorts are
%    stored. Spikes.man_sort is a vector of 0/1 for each channel, to be
%    used as a flag for which channels have been manually sorted. 
%
%    2 functions called:
%     - [thresh, reject] = calculate_thresholds(Phys)
%           Prompts user to select clean trials, for calculating thresholds
%           based on std for each channel.
%     - spks = pp_sort_channel
%           Currently, called individually for each channel. This script
%           is the one that calls the UMS functions to run the sorting. 
%    
%    Saves data files:
%     SUBJECT_sess-XX_Spikes.mat
%    
%  KP, 2016-04; last updated 2018-02
%  2018-02: Now written to work with TrialData table instead of SoundData
%    format for saving stimulus data. 
%

%%
close all
tic

global fs fn

% Set directories based on processor used 
fn = set_paths_directories(subject,session,1);

% Load data structures
fprintf('\nloading data...')
filename = sprintf('%s_sess-%s_Phys',subject,session);
load(fullfile(fn.processed,subject,filename));
filename = sprintf('%s_sess-%s_Info',subject,session);
load(fullfile(fn.processed,subject,filename));
filename = sprintf('%s_sess-%s_TrialData',subject,session);
load(fullfile(fn.processed,subject,filename));

fs = Info.fs;
fprintf(' done.\n')


%%

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Call program to automatically select clean data segments and set thresholds
[thresh, reject] = calculate_thresholds_auto( Phys, TrialData, Info);
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


%%
Spikes = struct();

for ich = 1:size(Phys,1)
    
    data = Phys(ich,:);
    fprintf('sorting ch %i... ',ich)
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    spks = pp_sort_channel(data, thresh(ich), reject(ich));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    Spikes.channel(ich)   = ich;
    Spikes.man_sort(ich)  = 0;
    Spikes.clustered(ich) = spks;
    
    clear spks
    
    % Save data after each channel
    savename = sprintf('%s_sess-%s_Spikes',subject,session);
    save(fullfile(fn.processed,subject,savename),'Spikes','-v7.3');
    
end

% Replicate the raw sorted data, where the final clus will be saved after
% manual sorting.
Spikes.sorted = Spikes.clustered;



% SAVE DATA 

% Save completed Spikes file
fprintf('\nsaving data...\n')
savename = sprintf('%s_sess-%s_Spikes',subject,session);
save(fullfile(fn.processed,subject,savename),'Spikes','-v7.3');



%
fprintf('\n**Finished sorting and saving data for this session.\n')
fprintf('------------------------------------------------------------\n')

toc
close all



end