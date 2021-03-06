function ProcessPhys_SynKS( SESS_LABEL )
%
%  pp_processPhys_Syn( subject, session_label )
%    
%    ** Trials with artifact not yet marked
%
%    Saves data files:
%     SUBJECT_sess-XX_Info.mat
%     SUBJECT_sess-XX_TrialData.mat
%     and a folder called "sorting", which contains kilosort files
%
%  KP, 2018-08
%

%%
        %     fs = epData.streams.rVrt.fs (1/5 of phys data)
        % epData.streams.rVrt.data
        %  (1,:) = Instantaneous AM rate  (of the period that is just ending)
        %  (2,:) = Sound output          
        %  (3,:) = AM depth
        %  (4,:) = dB SPL
        %  (5,:) = HP
        %  (6,:) = LP
        %  (7,:) = Spout TTL
        %  (8,:) = ITI TTL
        
        % falling edge of InTrial:
        %   epData.epocs.RCod  =  response code
        % StimTrial TTL:
        %   epData.epocs.TTyp  =  TrialType (0 safe, 1 warn)
        %   epData.epocs.Opto  =  optostim
        %   epData.epocs.rVID  =  rateVec_ID
        % rising edge of PHASE0:
        %   epData.epocs.AMrt  =  AMRATE
        
        
%%

global fn
close all
tic

% Set directories for loading and saving data
fn = set_paths_directories();


% Load this block epData datafile

[BLOCK, BLOCKPATH] = uigetfile(fn.raw, 'Select block file to process...');
foo = strsplit(BLOCKPATH,filesep);
SUBJECT = foo{end-2};

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
    keyboard  %need stimfiles
% Info.behav_fn      = epData.info.fnBeh;
Info.fs_sound      = epData.streams.rVrt.fs;
Info.sound_rows    = sound_rows;
% Info.stimfiles     = epData.stimfs;
end
% Info.noiseremvsamp = noisewin;   %for KS
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


%%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse stimulus data and save necessary info

if strcmp(BlockType,'behavior')

    [TrialData,SpoutStream,SoundStream,Info] = pp_parse_sound_data( epData.streams.rVrt.data, epData.epocs, Info );

end


%% SAVE DATA

fprintf(' saving Info struct and Trial data...')


saveDataDir = fullfile(fn.processed,SUBJECT);
if ~exist(saveDataDir,'dir')
    mkdir(saveDataDir)
end

% Save Info structure
savename = sprintf('%s_sess-%s_Info',SUBJECT,SESS_LABEL);
save(fullfile( saveDataDir, savename),'Info','-v7.3');


if strcmp(BlockType,'behavior') && exist('TrialData','var')
    
    % Save Stimulus data variables
    savename = sprintf('%s_sess-%s_TrialData',SUBJECT,SESS_LABEL);
    save(fullfile(saveDataDir, savename),'TrialData','SpoutStream','SoundStream','-v7.3');

end

fprintf('  done.\n')



%% Preprocess physiology data and run KiloSort algorithm

%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

fprintf(' removing noisy pre-session portion...')

% Remove beginning of session, before headstage & animal in room

data = epData.streams.RSn1.data(1,:);  %only need to use one channel -- if ch 1 is dead, change to another

begRMS = rms(data(1,1:round(Info.fs*10)));
endRMS = rms(data(1,(end-round(Info.fs*21)):(end-round(Info.fs*1)) ));

if (begRMS-1.2*endRMS)<0.01 %while testing
    keyboard
end

sampwin  = round(Info.fs*2); %2000 ms
sampstep = round(Info.fs/10); %500 ms

isamp = 1;
RMS_HIGH = 1;
while RMS_HIGH
    thisRMS = rms(data(1,isamp:(isamp+sampwin-1)));
    if thisRMS < 1.2*endRMS
        RMS_HIGH = 0;
    else
        isamp = isamp+sampstep;
    end
end

% % % Plot to check
% % begwins  = 1:sampstep:size(data,2);
% % roughRMS = nan(size(begwins));
% % for iw = 1:numel(begwins)
% %     roughRMS(iw) = rms(data(1,begwins(iw):min(begwins(iw)+sampwin-1,size(data,2))));
% % end
% % 
% % hf=figure;
% % plot(begwins,roughRMS,'k')
% % hold on
% % plot([isamp isamp],[0 begRMS],'r','LineWidth',2)
% % 
% % close(hf);

% Remove data in beginning
dataClean = epData.streams.RSn1.data;
dataClean(:,1:isamp) = 0;

fprintf('  done.\n')


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% Where to save the sorting files to
saveDataDir = fullfile(fn.sessdata,'sorting');
if ~exist(saveDataDir,'dir')
    mkdir(saveDataDir)
end

% 
% % Convert raw physiology data from epData to a .dat file
% fprintf(' writing dat file...')
datName = [SUBJECT '_' SESS_LABEL '_' strtok(BLOCK,'.') '.dat']; %must match string in 
fidn=fopen(fullfile(saveDataDir,datName),'w');
fwrite(fidn,int16( dataClean * 1e6),'int16'); %(int16 format requires integers)
fclose(fidn);  
fprintf('  done.\n')



%% SPIKE SORTING (kilosort)

fprintf('\n  ~~~~~~~~~ KILOSORT ~~~~~~~~~\n')


% Run configuration for this probe

% whichProbe = 'NN4x4';
% whichProbe = 'Atlas4tet';
whichProbe = 'Buzsaki5x12_64';

configFileDir = fullfile(fn.processed,'SortingConfig',whichProbe);
run(fullfile(configFileDir, ['config_' whichProbe '.m']))



% Run Kilosort processing and algorithm
[rez, DATA, uproj] = preprocessData_kp(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


% Save the data
rezToPhy_kp(rez);
delete(ops.fproc); %just deletes empty temp_wh.dat file


%%

%
fprintf('\n\n**Finished processing and saving data for this recording session.\n')
fprintf('------------------------------------------------------------\n\n')

%
toc

end




