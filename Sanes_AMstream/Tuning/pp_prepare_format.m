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
%    >> for AM stream expt, tuning blocks only
%    
%    For each block, calls functions:
%     - pp_make_stim_struct(epData)
%     - pp_make_phys_struct(epData,t1,t2)
%
%    Saves data files:
%     SUBJECT_sess-XX_Info.mat
%     SUBJECT_sess-XX_Phys.mat
%     SUBJECT_sess-XX_Stim.mat
%   
%  KP, 2016-04; last updated 2017-06
% 

%%

tic

% Set directories for loading and saving data
fn = set_paths_directories;


%% Process data...

fprintf('\n------------------------------------------------------------')
fprintf('\nProcessing ephys data from: %s',subject)



Stim=[];   Phys=[];   info_stimfiles = cell(1,numel(BLOCKS));

for ib=1:numel(BLOCKS)
    
    % Choose folder containing WAV files used in this block
    
    fprintf('\n BLOCK %s\n',num2str(BLOCKS(ib)))
    
        
    % Load this block datafile
    
    block_str = sprintf('Block-%i.mat',BLOCKS(ib));
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
    end
    starttimes{ib} = epData.info.starttime;
    
        
    
    % Create STIM struct with stimulus info
    try
        
    fprintf('\n   adding to STIM struct...')
    st=[]; epoc_trs=[];
    [st, epoc_trs] = pp_make_stim_struct( epData, BLOCKS(ib) );
    Stim = [Stim st];
    fprintf('                           done.')
    
    catch
        keyboard
    end
    
    
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % Set extraction window
    if ib==1
        
        t1ms = -299;         %ms of data to pull before tr onset
        t2ms =  400;         %ms of data to pull after tr onset
        
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


% Create Info struct
Info = struct;
Info.subject      = subject;
Info.date         = epData.info.date;
Info.starttimes   = starttimes;
Info.blocks       = BLOCKS';
Info.t_win_ms     = [t1ms t2ms];
Info.fs           = epData.streams.Wave.fs;
Info.artifact_trs = artifact_trs;



% SAVE DATA 

fprintf('\nsaving data...')

try
% Save Stim structure
savename = sprintf('%s_sess-%s_Stim',subject,session_label);
save(fullfile( fn.processed,subject,savename),'Stim','-v7.3');

% Save Phys structure
savename = sprintf('%s_sess-%s_Phys',subject,session_label);
save(fullfile( fn.processed,subject,savename),'Phys','-v7.3');

% Save Info structure
savename = sprintf('%s_sess-%s_Info',subject,session_label);
save(fullfile( fn.processed,subject,savename),'Info','-v7.3');

catch
    keyboard
end

%
fprintf('\n\n**Finished processing and saving data for this recording session.\n')
fprintf('------------------------------------------------------------\n\n')

%
toc

end  



