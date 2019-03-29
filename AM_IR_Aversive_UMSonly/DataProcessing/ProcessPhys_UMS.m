function ProcessPhys_UMS( BLOCK, subject, session_label, SKIP_PHYS )
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
%  KP, 2016-04; last updated 2018-01
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
tic

global fn

if numel(BLOCK)>1
    keyboard
end

close all

%% Process data...

% Set directories for loading and saving data
fn = set_paths_directories(subject,session_label,1);


fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from: %s',subject)


if ~isempty(BLOCK)
    
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
    
    
else %acquired with synapse
    
    keyboard
    
end


%%
% Determine recording type (aversive block or tuning block)

if isfield(epData.streams,'rVrt')
    BlockType = 'behavior';
else
    BlockType = 'tuning';
    keyboard
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
Info.subject       = subject;
Info.session       = session_label;
Info.date          = epData.info.date;
Info.block         = BLOCK';
Info.blockType     = BlockType;
Info.fs            = epData.streams.Wave.fs;
if strcmp(BlockType,'behavior')
Info.fs_sound      = epData.streams.rVrt.fs;
end
Info.sound_rows    = sound_rows;
Info.stimfiles     = epData.stimfs;
% Info.noiseremvsamp = noisewin;   %for KS
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


%%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse stimulus data and save necessary info

[TrialData,SpoutStream,SoundStream,Info] = pp_parse_sound_data( epData.streams.rVrt.data, epData.epocs, Info );


%%

if nargin<4 || ~SKIP_PHYS

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Filter PHYS data and keep in one long stream
try
    fprintf('\n   adding to PHYS struct...')
    Phys=[];
    Phys = pp_filter_phys_data( epData );
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
    savedir = fullfile(fn.processed,subject);
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    % Save Info structure
    savename = sprintf('%s_sess-%s_Info',subject,session_label);
    save(fullfile( savedir, savename),'Info','-v7.3');
    
    % Save Phys structure
    if nargin<4 || ~SKIP_PHYS
        savename = sprintf('%s_sess-%s_Phys',subject,session_label);
        save(fullfile( savedir, savename),'Phys','-v7.3');
    end
    
    % Save Stimulus data variables
    savename = sprintf('%s_sess-%s_TrialData',subject,session_label);
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

end




