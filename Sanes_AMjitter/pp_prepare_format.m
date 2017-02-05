
function pp_prepare_format( BLOCKS, subject, session_label )
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
%     - pp_make_stim_struct(epData)
%     - pp_make_phys_struct(epData,t1,t2)
%
%    Saves data files:
%     SUBJECT_sess-XX_Info.mat
%     SUBJECT_sess-XX_Phys.mat
%     SUBJECT_sess-XX_Stim.mat
%     SUBJECT_sess-XX_LFP.mat
%     SUBJECT_sess-XX_Wavlet.mat
%
%    Now working with data locally.
%   
%  KP, 2016-04; last updated 2016-10
% 

tic

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %To Add:
    %   change time window with stim duration?
    %   remove noisy trials altogether?
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
% Set directories for loading and saving data

includePath='/Users/kpenikis/Documents/MATLAB/ums2k_02_23_2012';
addpath(genpath(includePath));
addpath(genpath('/Users/kpenikis/Documents/MATLAB/eeglab13_6_5b/functions/'))
addpath('helpers')

% Raw data location
datadirlocal = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData';
datadirdrive = '/Volumes/Seagate Backup Plus Drive/RawData';


% Saving location
savedir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
if ~exist('session_label','var')
    prompt = '\n--> please enter uppercase letters for session label and press enter.';
    session_label = input(prompt,'s');
end

if ~exist(fullfile(savedir,subject),'dir')
    [~,~,~] = mkdir(fullfile(savedir,subject));
end


% Process data...

fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from: %s',subject)



Stim=[];   Phys=[];   info_stimfiles = cell(1,numel(BLOCKS));
data_lp3=[];  data_lp30=[]; data_lp400=[]; data_raw=[];

for ib=1:numel(BLOCKS)
    
    % Choose folder containing WAV files used in this block
    
    fprintf('\n BLOCK %s\n',num2str(BLOCKS(ib)))
    
        
    % Load this block datafile
    
    block_str = sprintf('Block-%i.mat',BLOCKS(ib));
    datafile = fullfile(datadirlocal,subject,block_str);
    if ~exist(datafile,'file') %if doesnt exist locally, check for harddrive
        datafile = fullfile(datadirdrive,subject,block_str);
        if ~exist(datafile,'file')
            warning('External hard drive not connected, and file not found locally.')
            keyboard
        end
    end
    fprintf(' loading data file %s...',datafile)
    clear epData;
    load(datafile,'-mat'); %loads data struct: epData
    if isempty(epData)
        error('data file did not load correctly!') 
    end
    starttimes{ib} = epData.info.starttime;
    
    
    % Get names of stimulus wavfiles
    
    if isfield(epData,'wavfilenames')
        stimfiles = epData.wavfilenames;
        
    elseif isfield(epData,'matfilenames')
%         stfns = dir(fullfile('/Users/kpenikis/Documents/SanesLab/Data/raw_data',subject,sprintf('Block-%i_Stim',BLOCKS(ib)),'*.mat'));
%         stfns = {stfns.name};
%         stimfiles = sort_fns(stfns);
%         clear stfns
        stimfiles = epData.matfilenames;
        
    else 
        stimfiles = 'basic_tuning';
    end
    
    info_stimfiles{ib} = stimfiles;
    
    
    % Create STIM struct with stimulus info
    try
        
    fprintf('\n   adding to STIM struct...')
    st=[]; epoc_trs=[];
    [st, epoc_trs] = pp_make_stim_struct( epData, BLOCKS(ib), stimfiles );
    Stim = [Stim st];
    fprintf('                           done.')
    
    catch
        keyboard
    end
    
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % Set extraction window
    if ib==1
        if strcmp(stimfiles,'basic_tuning')
            maxDur = 1150;
        else
            maxDur = max([Stim.stimDur]); %3300;  %soft coding only works if jitter block is first
        end
        
        t1ms = -799;         %ms of data to pull before tr onset
        t2ms = maxDur+500;   %ms of data to pull after tr onset
        
        fs = epData.streams.Wave.fs;
        
        t1 = round(t1ms/1000*fs);   %converted to samples
        t2 = round(t2ms/1000*fs);   %converted to samples
    end
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    % Create struct of PHYS data in windows around trials
    
    try
        
    fprintf('\n   adding to PHYS struct...')
    ph=[];
    ph = pp_make_phys_struct( epData, t1, t2 , epoc_trs);
    Phys = [Phys; ph];
    fprintf('       done.')
    
    catch
        keyboard
    end
    
% % %     % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% % %     
% % %     
% % %     % Create struct of LFP data in windows around trials
% % %     
% % %     try
% % %         fprintf('\n   adding to LFP struct... ')
% % %         [ lp3, lp30, lp400, raw, fs_lfp ] = pp_make_lfp_struct(epData,t1,t2); 
% % %         data_lp3  = [data_lp3; lp3];
% % %         data_lp30  = [data_lp30; lp30];
% % %         data_lp400 = [data_lp400; lp400];
% % %         data_raw   = [data_raw; raw];
% % %         lp3=[]; lp30=[]; lp400=[]; raw=[];
% % %         
% % %     catch
% % %         keyboard
% % %     end
% % %         
% % %     
% % %     % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
% % %     
% % %     % Create struct of LFP WAVLET data in windows around trials
% % % 
% % %     t1 = round(t1ms/1000*fs);   %converted to samples
% % %     t2 = round(t2ms/1000*fs);   %converted to samples
% % %     
% % %     if numel(BLOCKS)>1
% % %         disp('Skipping wavlet analysis. Cant yet handle combining data across blocks.')
% % %         continue
% % %     else
% % %     
% % %     try
% % %         fprintf('\n   adding to Wavelet struct... ')
% % %         
% % %         for iband = [1 26 51 76]
% % %             
% % %             clear Wavlet
% % %             [Wavelet,fs_wavlet] = pp_make_wavlet_struct(epData, t1ms, t2ms, iband:iband+24);
% % %             
% % %             % Save Wavlet structure
% % %             wavlet_savename = sprintf('%s_sess-%s_Wavlets%i-%i',subject,session_label,iband,iband+24);
% % %             save(fullfile(savedir,subject,wavlet_savename),'Wavelet','-v7.3');
% % %             
% % %         end
% % %     catch
% % %         keyboard
% % %     end
% % %     
% % %     end
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
end  %for ib=1:numel(blocks)


% Find trials with artifact
[artifact_trs] = get_artifacts(Phys, fs);


% test_downsampling(data_lp400, fs_lfp, [t1ms t2ms], artifact_trs)

% % % % Combine lfp data into struct
% % % LFPdata.lp3    = data_lp3;
% % % LFPdata.lp30   = data_lp30;
% % % LFPdata.lp400  = data_lp400;
% % % LFPdata.fs     = fs_lfp;



% Create Info struct
Info = struct;
Info.subject      = subject;
Info.date         = epData.info.date;
Info.starttimes   = starttimes;
Info.blocks       = BLOCKS';
Info.t_win_ms     = [t1ms t2ms];
Info.fs           = epData.streams.Wave.fs;
% Info.fs_lfp       = fs_lfp;
Info.artifact_trs = artifact_trs;
Info.stimfiles    = info_stimfiles;



% SAVE DATA 

fprintf('\nsaving data...')

try
% Save Stim structure
savename = sprintf('%s_sess-%s_Stim',subject,session_label);
save(fullfile(savedir,subject,savename),'Stim','-v7.3');

% Save Phys structure
savename = sprintf('%s_sess-%s_Phys',subject,session_label);
save(fullfile(savedir,subject,savename),'Phys','-v7.3');

% Save LFP structure
% savename = sprintf('%s_sess-%s_LFP',subject,session_label);
% save(fullfile(savedir,subject,savename),'LFPdata','-v7.3');

% Save Info structure
savename = sprintf('%s_sess-%s_Info',subject,session_label);
save(fullfile(savedir,subject,savename),'Info','-v7.3');

catch
    keyboard
end

%
fprintf('\n\n**Finished processing and saving data for this recording session.\n')
fprintf('------------------------------------------------------------\n\n')

%
toc

end  

function sorted_fns = sort_fns(fns)

[~,idx] = sort(str2double(strtok(extractBetween(fns,'4Hz_','.mat'),'_')));
sorted_fns = {fns{idx}};

end



