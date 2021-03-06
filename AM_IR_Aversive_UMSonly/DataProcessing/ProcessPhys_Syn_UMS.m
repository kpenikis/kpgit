function ProcessPhys_Syn_UMS( SESS_LABEL, SKIP_PHYS )
%
%  pp_prepare_format( BLOCKS, subject, session_label )
%    Creates data structures in the desired format for processing. Groups
%    together sequentially-recorded blocks, indicated in input argument,
%    and makes a structure with stimulus info by trial, a structure with
%    ephys data in windows around each trial (in the format needed for
%    sorting in UMS2000) and file with relevant information, including the
%    time window around trials start for which ephys data was extracted.
%
%    Saves data files:
%     SUBJECT_sess-XX_Info.mat
%     SUBJECT_sess-XX_Phys.mat
%     SUBJECT_sess-XX_TrialData.mat
%
%  KP, 2018-12
%


%%
        % fs = epData.streams.rVrt.fs 
        % epData.streams.rVrt.data
        %
        % falling edge of InTrial:
        %   epData.epocs.RCod  =  response code
        %
        % StimTrial TTL:
        %   epData.epocs.TTyp  =  TrialType (0 safe, 1 warn)
        %   epData.epocs.Opto  =  optostim
        %   epData.epocs.rVID  =  rateVec_ID
        %
        % rising edge of PHASE0:
        %   epData.epocs.AMrt  =  AMRATE
        %

        
%%

global fn
close all
tic

% Set directories for loading and saving data
fn = set_paths_directories();


% Load this block epData datafile

[BLOCK, BLOCKPATH] = uigetfile(fn.raw, 'Select block file to process...');
foo = strsplit(BLOCKPATH,filesep);
SUBJECT = foo{end-1};

fn = set_paths_directories(SUBJECT,SESS_LABEL);


fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from BLOCK: %s\n',BLOCK)


%% Load datafile...

datafile = fullfile(BLOCKPATH,BLOCK);      %change back to fn.raw
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

BLOCK = strtok(BLOCK,'.');


%%
% Determine recording type (aversive block or tuning block)

if isfield(epData.streams,'rVrt')
    BlockType = 'behavior';
else
    BlockType = 'tuning';
end


%%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Create INFO struct

%%% CHECK FOR EARLY OR LATE RECORDING
if isfield(epData.epocs,'LPxx')
    sound_rows = {'Instantaneous AM rate';...
        'Sound output';...
        'AM depth';...
        'dB SPL';...
        'Phase0';...
        'StimTrial';...
        'Spout TTL';...
        'ITI' };
else
    sound_rows = {'Instantaneous AM rate';...
        'Sound output';...
        'AM depth';...
        'dB SPL';...
        'HP';...
        'LP';...
        'Spout TTL';...
        'ITI' };
end

Info = struct;
Info.subject       = SUBJECT;
Info.session       = SESS_LABEL;
Info.date          = epData.info.date;
Info.blockType     = BlockType;
Info.epData_fn     = BLOCK;
Info.fs            = epData.streams.RSn1.fs;
if strcmp(BlockType,'behavior')
Info.behav_fn      = epData.info.fnBeh;
Info.fs_sound      = epData.streams.rVrt.fs;
Info.sound_rows    = sound_rows;
Info.stimfiles     = epData.stimfs;
end

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


%%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse stimulus data and save necessary info

[TrialData,SpoutStream,SoundStream,Info] = pp_parse_sound_data( epData.streams.rVrt.data, epData.epocs, Info );


%%

if nargin<2 || ~SKIP_PHYS

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Process PHYS data
try
    fprintf('\n   processing PHYS data...')
    
    % Check if extended noisy period in beginning to exclude
    SampBegin = check_RMS( epData.streams.RSn1.data(8,:), epData.streams.RSn1.fs );
    if SampBegin>1
        epData.streams.RSn1.data(:,1:SampBegin) = 0;
    end
    
    % Filter and whiten data
    Phys=[];
    Phys = pp_filter_phys_data( double(epData.streams.RSn1.data), epData.streams.RSn1.fs );
    
    % Synapse saves raw raw data and UMS does not whiten, so now apply
    % whitening or common average referencing. 
    useChs = [1:7 9:16];
    Phys = pp_common_avg_ref_data( Phys, useChs );
    
    fprintf('       done.\n')
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Find trials with artifacts
 Info = identify_artifact_trials(Info,Phys,TrialData);
 % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

catch
    keyboard
end

end %SKIP_PHYS



%% SAVE DATA

fprintf('\nsaving data...')

try
    savedir = fullfile(fn.processed,SUBJECT);
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    % Save Info structure
    savename = sprintf('%s_sess-%s_Info',SUBJECT,SESS_LABEL);
    save(fullfile( savedir, savename),'Info','-v7.3');
    
    % Save Phys structure
    if nargin<2 || ~SKIP_PHYS
        savename = sprintf('%s_sess-%s_Phys',SUBJECT,SESS_LABEL);
        save(fullfile( savedir, savename),'Phys','-v7.3');
    end
    
    % Save Stimulus data variables
    savename = sprintf('%s_sess-%s_TrialData',SUBJECT,SESS_LABEL);
    save(fullfile(savedir, savename),'TrialData','SpoutStream','SoundStream','-v7.3');
    
    %
    
catch
    keyboard
end

%
fprintf('\n\n**Finished processing and saving data for this recording session.\n')
fprintf('------------------------------------------------------------\n\n')

%
toc


%%

pp_sort_session( SUBJECT, SESS_LABEL, Phys, Info, TrialData )


end




