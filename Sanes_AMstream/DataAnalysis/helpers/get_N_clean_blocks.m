function [N,minTime] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd)
% Get number of clean trials for each block type, with these sound
% parameters. Calls function "get_blockOnsets". 
%  KP, 207-08

N = nan(1,10);
time = nan(1,10);

% Go through each stimulus
for ib = 1:numel(Info.blockKey)
    
    %if unmodulated or silent, skip for now
    if ib>=11, continue, end
    
    [~,~,bkStart_ms,bkStop_ms] =  get_blockOnsets( SoundData,ib,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
    
    N(ib) = numel(bkStart_ms);
    
    time(ib) = min(bkStop_ms-bkStart_ms);
    
end

minTime = min(time);

end