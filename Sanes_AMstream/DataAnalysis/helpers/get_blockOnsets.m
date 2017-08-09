function [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms] = get_blockOnsets(SoundData,...
                                            ib,spl,lpn,amd,ArtifactFlag,fs)
% [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms] = get_blockOnsets(SoundData,...
%                                            ib,spl,lpn,amd,ArtifactFlag,fs)
%   Gets the start and stop times (in samples and ms) for a given block
%   with a given set of sound parameters. 
%   Removes all blocks that are not the specified spl, lp noise, am depth,
%   OR that contain artifact 
%   OR that animal was not on the spout
%   OR are too short in duration (when recording cut off early)
%
% KP, 2017-08
%

                                        
% Get samples of beginning and end of block
bkStart_samps = 1+find( diff(SoundData(8,:)==ib) ==  1 );
bkStop_samps  =   find( diff(SoundData(8,:)==ib) == -1 );

% Remove all blocks that are not this spl, lp noise, am depth,
% OR that contain artifact 
% OR that animal wasnt on the spout
% OR are too short in duration (when recording cut off early)
rm_bk = [];
for it = 1:numel(bkStart_samps)
    if ib==5 && (SoundData(8,bkStart_samps(it)-1)==11)
        rm_bk = [rm_bk it];
    elseif ib~=5 && (SoundData(8,bkStart_samps(it)-1)==11)
        keyboard
    end
    if any(intersect(bkStart_samps(it):bkStop_samps(it),ArtifactFlag))...
            || ~all(SoundData(3, bkStart_samps(it):bkStop_samps(it)) ==amd)...
            || ~all(SoundData(4, bkStart_samps(it):bkStop_samps(it)) ==spl)...
            || ~all(SoundData(6, bkStart_samps(it):bkStop_samps(it)) ==lpn)...
            || ~all(SoundData(7, bkStart_samps(it):bkStop_samps(it)) ==1)...
            || ((bkStop_samps(it)-bkStart_samps(it))/fs) < 1.9
        rm_bk = [rm_bk it];
    end
end
bkStart_samps(rm_bk) = [];
bkStop_samps(rm_bk) = [];

% Convert to ms
bkStart_ms = round( bkStart_samps / fs*1000 );
bkStop_ms  = round( bkStop_samps  / fs*1000 );


if numel(bkStart_ms) ~= numel(bkStop_ms)
    keyboard
end


end