function [N,minTime] = get_N_clean_blocks_SS(SoundData,Info,ArtifactFlag,spl,amd,these_blocks)
% Get number of clean trials for each block type, with these sound
% parameters. Calls function "get_blockOnsets". 
%  KP, 207-08

N = nan(1,10);
time = nan(1,10);

if nargin<7 && ~exist('these_blocks','var')
    these_blocks = 1:10;
end

% Go through each stimulus
for ib = these_blocks
    
    %if unmodulated or silent, skip for now
    if ib>=11, continue, end
    
    [~,~,bkStart_ms,bkStop_ms] =  get_blockOnsets_SS( SoundData,ib,spl,amd,ArtifactFlag,Info.fs_sound);
    
    N(ib) = numel(bkStart_ms);
    
    if ~isempty(bkStart_ms)
        time(ib) = min(bkStop_ms-bkStart_ms);
    end
    
end

minTime = min(time);

end