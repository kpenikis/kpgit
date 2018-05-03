function [SoundData,Info] = pp_parse_sound_stream_v2(SoundData,Info)
%
% called by: pp_make_InfoPhys
% 


global fn

stimfolder = strsplit(Info.stimfiles, '\');
stimfolder = stimfolder{end};
stimfolder = fullfile(fn.stim,stimfolder);

stimfiles = dir(fullfile(stimfolder,'*.mat'));

if numel(stimfiles)==2
    protocol = 'Linearity';
elseif numel(stimfiles)==3
    protocol = 'SpectralSwitch';
elseif numel(stimfiles)==311
    protocol = 'Trials';
else
    keyboard
end

switch protocol
    
    case 'Linearity'
        
        % Load in the relevant stimulus files
        q=load(fullfile(stimfolder,'fullStream_AMrateVec_26trs_20160323.mat'));
        pd_AMrates = q.buffer;
        clear q
        
        q=load(fullfile(stimfolder,'fullStream_BlockVec_26trs_20160323.mat'));
        pd_Blocks  = q.buffer;
        blockKey   = q.stim_array;
        clear q
        
        % Call program to assign block labels
        [SoundData,Info] = get_block_labels(pd_AMrates,pd_Blocks,blockKey,SoundData,Info);
        
        
    case 'SpectralSwitch'
        
        % Load in the relevant stimulus files
        q=load(fullfile(stimfolder,'fullStream_AMrateVec_abriged.mat'));
        pd_AMrates = q.buffer;
        clear q
        
        q=load(fullfile(stimfolder,'fullStream_BlockVec_abriged.mat'));
        pd_Blocks  = q.buffer;
        blockKey   = q.stim_array;
        clear q
        
        q=load(fullfile(stimfolder,'fullStream_TargetVec_abriged'));
        currentCF  = q.buffer;
        clear q
        
        % Add an empty row
        SoundData = [SoundData; nan(1,size(SoundData,2))];
        
        % Call program to assign block labels
        [SoundData,Info] = get_block_labels(pd_AMrates,pd_Blocks,blockKey,SoundData,Info);
        
        
            
        
    case 'Trials'
        keyboard
        
end




end






