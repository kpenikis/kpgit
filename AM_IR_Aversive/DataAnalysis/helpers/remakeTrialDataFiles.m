function remakeTrialDataFiles(SUBJECT,SESSION)
% Differences: 
%   ITI spout pct calculated just in first 1.5 s now. 
%   New variable: RateStream
% 

fn = set_paths_directories(SUBJECT,SESSION);

% Load Info file
filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));


% Select block of raw data

[BLOCK, BLOCKPATH] = uigetfile(fn.raw, 'Select block file to process...');
foo = strsplit(BLOCKPATH,filesep);
SUBJECT = foo{end-1};

saveDir = fullfile(fn.processed,SUBJECT);
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from BLOCK: %s\n',BLOCK)
fprintf(' SUBJECT:    %s\n',SUBJECT)
fprintf(' SESS_LABEL: %s\n',SESSION)
fprintf(' save location: %s\n',saveDir)


% Load epData file

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


%% Call parse sound data function

[TrialData,SpoutStream,SoundStream,RateStream,Info] = pp_parse_sound_data( epData.streams.rVrt.data, epData.epocs, Info );


%% SAVE DATA

savename = sprintf('%s_sess-%s_TrialData',SUBJECT,SESSION);
save(fullfile(saveDir, savename),'TrialData','SpoutStream','SoundStream','RateStream','-v7.3');

fprintf('New session data saved.\n')



end