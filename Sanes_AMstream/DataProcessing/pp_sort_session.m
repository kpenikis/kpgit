function Spikes = pp_sort_session( subject, session, Spikes )
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
%    pp_sort_session( ..., Spikes )  load a Spikes struct into the
%    workspace and include as an input argument in order to begin sorting
%    from a later channel. 
%
%  KP, 2016-04; last updated 2017-04
%

%%
close all
tic

global fs fn

% Set directories based on processor used 
fn = set_paths_directories;

% Load data structures
fprintf('\nloading data...')
filename = sprintf('%s_sess-%s_Phys',subject,session);
load(fullfile(fn.processed,subject,filename));
filename = sprintf('%s_sess-%s_Info',subject,session);
load(fullfile(fn.processed,subject,filename));

if isfield('Info','t_win_ms') %diff(Info.t_win_ms)<1000
    filename = sprintf('%s_sess-%s_Stim',subject,session);
    load(fullfile(fn.processed,subject,filename));
else
    filename = sprintf('%s_sess-%s_SoundData',subject,session);
    load(fullfile(fn.processed,subject,filename));
end

fs = Info.fs;
fprintf(' done.\n')


%%

if exist('SoundData','var')
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % Call program to automatically select clean data segments and set thresholds
    [thresh, reject] = calculate_thresholds_auto(Phys, SoundData, Info);
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

else
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % Manually select which flagged trials contain disruptive artifact
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if ~any(strcmp(fieldnames(Info.artifact_trs),'manual'))
        
        Info = manually_mark_trials(Info,session,Stim,Phys);
        
        % Re-save Info structure
        filename = sprintf( '%s_sess-%s_Info',subject,session);
        save(fullfile(fn.processed,subject,filename),'Info','-v7.3');
        
        clear Stim
    end
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % And calculate thresholds
    
    segment_length_s = diff(Info.t_win_ms)/1000;  %length of each data segment (sec)
    [thresh, reject] = calculate_thresholds(Phys, segment_length_s);

    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

end

%%
Spikes = struct();

for ich = 1:size(Phys,1)
    if ich==8
        Spikes.channel(ich) = ich;
        continue
    end
    
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