
function fixStimInfoFilenameSorting( BLOCKS, subject, session_label )
%   
%  KP, 2016-11
% 

    
% Set directories for loading and saving data

includePath='/Users/kpenikis/Documents/MATLAB/ums2k_02_23_2012';
addpath(genpath(includePath));
addpath(genpath('/Users/kpenikis/Documents/MATLAB/eeglab13_6_5b/functions/'))


savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
if ~exist('session_label','var')
    prompt = '\n--> please enter uppercase letters for session label and press enter.';
    session_label = input(prompt,'s');
end

if ~exist(fullfile(savedir,subject),'dir')
    [~,~,~] = mkdir(fullfile(savedir,subject));
end


fprintf('\n------------------------------------------------------------')
fprintf('\nGetting ephys data from: %s',subject)

Stim=[];   info_stimfiles = cell(1,numel(BLOCKS));

for ib=1:numel(BLOCKS)
    
    % Choose folder containing WAV files used in this block
    
    fprintf('\n BLOCK %s\n',num2str(BLOCKS(ib)))
    
        
    % Load this block datafile
    
    block_str = sprintf('Block-%i.mat',BLOCKS(ib));
    datafile = fullfile('/Users/kpenikis/Documents/SanesLab/Data/raw_data',subject,block_str);
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
        stfns = dir(fullfile('/Users/kpenikis/Documents/SanesLab/Data/raw_data',subject,sprintf('Block-%i_Stim',BLOCKS(ib)),'*.mat'));
        stfns = {stfns.name};
        stimfiles = sort_fns(stfns);
        
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
        
    
end  %for ib=1:numel(blocks)

% Load Info struct and replace stimfile field
infoname = sprintf('%s_sess-%s_Info',subject,session_label);
load(fullfile(savedir,subject,infoname))

Info.stimfiles    = info_stimfiles;


% SAVE DATA 

fprintf('\nsaving data...')

try
    session_label = 'EA2';
% Save Stim structure
savename = sprintf('%s_sess-%s_Stim',subject,session_label);
save(fullfile(savedir,subject,savename),'Stim','-v7.3');

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
keyboard
end



