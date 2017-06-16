function pp_make_InfoPhys( BLOCK, subject, session_label )
%
%  pp_prepare_format( BLOCKS, subject, session_label )
%    Creates data structures in the desired format for processing. Groups
%    together sequentially-recorded blocks, indicated in input argument,
%    and makes a structure with stimulus info by trial, a structure with
%    ephys data in windows around each trial (in the format needed for
%    sorting in UMS2000) and file with relevant information, including the
%    time window around trials start for which ephys data was extracted.
%
%    For each block, calls functions:
%     - pp_make_phys_struct(epData,t1,t2)
%
%    Saves data files:
%     SUBJECT_sess-XX_Info.mat
%     SUBJECT_sess-XX_Phys.mat
%     SUBJECT_sess-XX_SoundData.mat
%
%  KP, 2016-04; last updated 2017-06
%

%%
        %     fs = epData.streams.rVrt.fs (1/5 of phys data)
        % epData.streams.rVrt.data
        %  (1,:) = Instantaneous AM rate <-- if Trials stim set, just this
        %  (2,:) = Sound output          <-- if Trials stim set, just this
        %  (3,:) = AM depth
        %  (4,:) = dB SPL
        %  (5,:) = HP
        %  (6,:) = LP
        %  (7,:) = Spout TTL
        
%%
tic

global fn

if numel(BLOCK)>1
    keyboard
end

close all

%% Process data...

% Set directories for loading and saving data
fn = set_paths_directories;


fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from: %s',subject)


% Choose folder containing WAV files used in this block

fprintf('\n BLOCK %s\n',num2str(BLOCK))


% Load this block epData datafile
block_str = sprintf('Block-%i.mat',BLOCK);
datafile = fullfile(fn.raw,subject,block_str);
if ~exist(datafile,'file')
    warning('raw epData file not found.')
    keyboard
end
fprintf(' loading data file %s...',datafile)
clear epData;
load(datafile,'-mat'); %loads data struct: epData
if isempty(epData)
    error('data file did not load correctly!')
else
    disp('done.')
end


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Create INFO struct
Info = struct;
Info.subject      = subject;
Info.session      = session_label;
Info.date         = epData.info.date;
Info.block        = BLOCK';
Info.fs           = epData.streams.Wave.fs;
Info.fs_sound     = epData.streams.rVrt.fs;
Info.stimfiles    = epData.info.stimpath;
Info.sound_rows   = {'Instantaneous AM rate';...
                     'Sound output';...
                     'AM depth';...
                     'dB SPL';...
                     'HP';...
                     'LP';...
                     'Spout TTL' };
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


if ~strcmp(epData.info.stimpath,'IR_AM_trials')

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Filter PHYS data and keep in one long stream
try
    fprintf('\n   adding to PHYS struct...')
    Phys=[];
    Phys = pp_filter_phys_data( epData );
    fprintf('       done.')
catch
    keyboard
end
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse and save SOUNDDATA matrix

SoundData = epData.streams.rVrt.data;

[SoundData,Info] = pp_parse_sound_stream_v2(SoundData,Info);
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Find windows with artifacts
[Info,SoundData] = identify_artifact_regions(Info,Phys,SoundData);
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Label clean blocks and calculate threshold for sorting
% [Info.thresh, Info.reject] = calculate_thresholds(Phys, SoundData, Info);
%%% --> now this will be done automatically at the start of spike sorting
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


else %if trials
    keyboard

end


%% SAVE DATA

fprintf('\nsaving data...')

try
    savedir = fullfile(fn.processed,subject);
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    % Save Info structure
    savename = sprintf('%s_sess-%s_Info',subject,session_label);
    save(fullfile( savedir, savename),'Info','-v7.3');
    
    % Save Phys structure
    savename = sprintf('%s_sess-%s_Phys',subject,session_label);
    save(fullfile( savedir, savename),'Phys','-v7.3');
    
    % Save SoundData structure
    savename = sprintf('%s_sess-%s_SoundData',subject,session_label);
    save(fullfile( savedir, savename),'SoundData','-v7.3');
    
    %
    
catch
    keyboard
end

%
fprintf('\n\n**Finished processing and saving data for this recording session.\n')
fprintf('------------------------------------------------------------\n\n')

%
toc

end




