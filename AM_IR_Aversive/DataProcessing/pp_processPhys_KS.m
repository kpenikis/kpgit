function pp_processPhys( BLOCK, subject, session_label )
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% IN PROGRESS, ADDING BLOCK IDS TO SD, SWITCHED TO UMS PROGRAM TO EDIT
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
%  KP, 2016-04; last updated 2018-01
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
tic

global fn

if numel(BLOCK)>1
    keyboard
end

close all

%% Process data...

% Set directories for loading and saving data
fn = set_paths_directories(subject,session_label);


fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from: %s',subject)

fprintf('\n BLOCK %s\n',num2str(BLOCK))


% Load this block epData datafile
block_str = sprintf('Block-%i.mat',BLOCK);
datafile = fullfile(fn.raw,subject,block_str);      %change back to fn.raw
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

sound_rows = {'Instantaneous AM rate';...
    'Sound output';...
    'AM depth';...
    'dB SPL';...
    'HP';...
    'LP';...
    'Spout TTL';...
    'ITI' };

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
% Info.noiseremvsamp = noisewin;   %dont remember what this was for
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


%%
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse and save SOUNDDATA matrix

SoundData = epData.streams.rVrt.data;

% onset of trials
tr_onset_rv  = epData.epocs.rVID.onset'; %MAKE SURE StimTrial GOES HIGH AT PHASE 0 OF FIRST AM PERIOD
% offset of trials
tr_offset_sd = (find(diff(SoundData(8,:))==1)/Info.fs_sound);
% stim ids
tr_id_ep = epData.epocs.rVID.data';

% If recording was stopped in the middle of a trial these sizes may not match
if size(tr_onset_rv,2) ~= size(tr_onset_rv,2) 
    keyboard
elseif any((tr_offset_sd-tr_onset_rv)<0) || any((tr_offset_sd-tr_onset_rv)>3)
    keyboard
end


% ITI rates
ITI_rates = unique(SoundData(1,2+find(diff(SoundData(8,:))==1)));
ITI_blockIDs = size(Info.stimfiles,1) + [1 2];
for ir = 1:numel(ITI_rates)
    Info.stimfiles{end+1} = ['ITI_' num2str(ITI_rates(ir))];
end

% % % SECONDS
% % figure; 
% % plot((1:1:size(SoundData,2))/Info.fs_sound,SoundData(8,1:1:end),'r') %iti
% % hold on
% % plot((1:1:size(SoundData,2))/Info.fs_sound,SoundData(2,1:1:end)/2+0.5,'Color',[0.6 0.6 0.6]) %sound
% % plot([tr_onset_rv' tr_offset_sd']',[tr_id_ep' tr_id_ep']'./max(tr_id_ep),'-k') %trial rvids
% % plot(epData.epocs.AMrt.onset,epData.epocs.AMrt.data./max(epData.epocs.AMrt.data),'g') %AM rate (stored at end of period)
% % plot((1:1:size(SoundData,2))/Info.fs_sound, SoundData(1,1:1:end)./max(SoundData(1,1:1:end)),'c') %AM rate (sound data)
% % title('SECONDS')
% % 
% % % SAMPLES
% % figure; 
% % plot((1:1:size(SoundData,2)),SoundData(8,1:1:end),'r') %iti
% % hold on
% % plot((1:1:size(SoundData,2)),SoundData(2,1:1:end)/2+0.5,'Color',[0.6 0.6 0.6]) %sound
% % plot([tr_onset_rv' tr_offset_sd']'.*Info.fs_sound, [tr_id_ep' tr_id_ep']'./max(tr_id_ep),'-k') %trial rvids
% % plot(epData.epocs.AMrt.onset.*Info.fs_sound, epData.epocs.AMrt.data./max(epData.epocs.AMrt.data),'g') %AM rate (stored at end of period)
% % plot((1:1:size(SoundData,2)), SoundData(1,1:1:end)./max(SoundData(1,1:1:end)),'c') %AM rate (sound data)
% % title('SAMPLES')
% % 


% End goal: add row to SoundData with trial code. rvID plus ITIs, plus 0
% for silence, plus NaN when not an analysis period (e.g. sound cut out, 
% extra period after iti but before next block began)

% Empty data vector
SD_addRow = nan(1,size(SoundData,2));

% Label trial blocks
for ib = 1:numel(tr_id_ep)
    ib_samps = ceil([tr_onset_rv(ib) tr_offset_sd(ib)] .*Info.fs_sound);
    SD_addRow(ib_samps(1):ib_samps(end)) = tr_id_ep(ib) * ones(1, diff(ib_samps)+1);
end

% Label ITI periods
if ~all(ITI_rates == unique(SoundData(1,2+find(diff(SoundData(8,:))==1))))
    keyboard
end

each_ITI_start = 1+find(diff(SoundData(8,:))==1);
each_ITI_end   = find(diff(SoundData(8,:))==-1);
if numel(each_ITI_start) ~= numel(each_ITI_end), keyboard, end

for ib = 1:numel(each_ITI_start)
    this_ITI = ITI_blockIDs( ITI_rates == double(SoundData(1,1+each_ITI_start(ib))) );
    SD_addRow(each_ITI_start(ib):each_ITI_end(ib)) = this_ITI * ones(1,numel(each_ITI_start(ib):each_ITI_end(ib)));
end
% Samples between end of ITI and start of next block
% round(tr_onset_rv(2:end).*Info.fs_sound - each_ITI_end(1:end-1))

% Label silence (when also drinking)
SD_addRow(1, SoundData(2,:)==0 & SoundData(7,:)==1 ) = 0;

figure(20);
plot(SD_addRow,'b')

SoundData = [SoundData; SD_addRow];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% IN PROGRESS, AND SWITCHED TO UMS PROGRAM TO EDIT
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Outstanding issues 
%  1) how to easily identify periodic-periodic transitions (iti to rvid
%  when there is a discontinuity/NaN separating them)
%  2) check that silent period is at the beginning
%  3) what to do if the number of ttls don't line up



% [SoundData,Info] = pp_parse_sound_stream_v2(SoundData,Info);

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

keyboard


%% Preprocess physiology data

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Filter data 
% dataFilt = pp_filter_phys_data(epData);

% Clean up artifacts
[Info,dataClean,SoundData] = identify_artifact_regions(Info,epData.streams.Wave.data,SoundData);
% [Info,dataClean,SoundData] = identify_artifact_regions(Info,dataFilt,SoundData);

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~



%% SPIKE SORTING (kilosort)

% Where to save the sorting files to
saveDataDir = fullfile(fn.sessdata,'sorting');
if ~exist(saveDataDir,'dir')
    mkdir(saveDataDir)
end

% Convert raw physiology data from epData to a .dat file
datName = [subject '_' session_label '_block' num2str(BLOCK) '_clean.dat'];
fidn=fopen(fullfile(saveDataDir,datName),'w');
fwrite(fidn,int16( dataClean * 1e6),'int16'); %(int16 format requires integers)
fclose(fidn);


% Check for probe geometry and config files
whichProbe = 'NN4x4';
whichProbe = 'Atlas4tet';
switch whichProbe
    case 'NN4x4'
        configFileDir = fullfile(fn.processed,'SortingConfig','NN4x4');
    case 'Atlas4tet'
        configFileDir = fullfile(fn.processed,'SortingConfig','Atlas4tet');
end

run(fullfile(configFileDir, ['config_' whichProbe 'probe.m']))


% Run Kilosort processing and algorithm
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


% Save the data
rezToPhy_kp(rez);
delete(ops.fproc);


 
return






%%

% Determine recording type
if isfield(epData.streams,'rVrt') && isfield(epData,'stimfs')
    stimfiles = epData.stimfs;
    datafn = 'rVrt';
    sound_rows = {'Instantaneous AM rate';...
                     'Sound output';...
                     'AM depth';...
                     'dB SPL';...
                     'HP';...
                     'LP';...
                     'Spout TTL';...
                     'ITI' };
% elseif isfield(epData.streams,'rVcf')
%     datafn = 'rVcf';
%     sound_rows = {'Instantaneous AM rate';...
%                      'Sound output';...
%                      'AM depth';...
%                      'dB SPL';...
%                      'CF';...
%                      'Spout TTL';...
%                      'placeholder'};
else
    warning('Might be processing the wrong data file')
    keyboard
end

% 
% % Determine stimulus set
% switch epData.info.stimset
%     case 'AC'
%         stimfiles = {'rateVec_2_AC.mat' 'rateVec_4_AC.mat' 'rateVec_8_AC.mat' 'rateVec_16_AC.mat' 'rateVec_32_AC.mat' };
%     case 'DB'
%         stimfiles = {'rateVec_2_DB.mat' 'rateVec_4_DB.mat' 'rateVec_8_DB.mat' 'rateVec_16_DB.mat' 'rateVec_32_DB.mat' };
%     case 'ACDB'
%         stimfiles = {'rateVec_2_AC.mat' 'rateVec_4_AC.mat' 'rateVec_8_AC.mat' 'rateVec_16_AC.mat' 'rateVec_32_AC.mat' 'rateVec_2_DB.mat' 'rateVec_4_DB.mat' 'rateVec_8_DB.mat' 'rateVec_16_DB.mat' 'rateVec_32_DB.mat' };
% 
% end



% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Create INFO struct
Info = struct;
Info.subject       = subject;
Info.session       = session_label;
Info.date          = epData.info.date;
Info.block         = BLOCK';
Info.fs            = epData.streams.Wave.fs;
Info.fs_sound      = epData.streams.(datafn).fs;
Info.sound_rows    = sound_rows;
Info.stimfiles     = epData.stimfs;
Info.noiseremvsamp = noisewin;
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~



try

% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Parse and save SOUNDDATA matrix

SoundData = epData.streams.(datafn).data;

[SoundData,Info] = pp_parse_sound_stream_v2(SoundData,Info);
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Find windows with artifacts
[Info,SoundData] = identify_artifact_regions(Info,Phys,SoundData);
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

catch
    keyboard
end
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% Label clean blocks and calculate threshold for sorting
% [Info.thresh, Info.reject] = calculate_thresholds(Phys, SoundData, Info);
%%% --> now this will be done automatically at the start of spike sorting
% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


keyboard

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




