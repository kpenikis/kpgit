function [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms] = get_blockOnsets_SS(SoundData,...
                                            ib,spl,amd,ArtifactFlag,fs)
% [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms] = get_blockOnsets_SS(SoundData,...
%                                            ib,spl,amd,ArtifactFlag,fs)
%   Gets the start and stop times (in samples and ms) for a given block
%   with a given set of sound parameters. 
%   Removes all blocks that are not the specified spl, lp noise, am depth,
%   OR that contain artifact 
%   OR that animal was not on the spout
%   OR are too short in duration (when recording cut off early)
%
% KP, 2017-08
%

scrsz = get(0,'ScreenSize');
figsize = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];
                                        
% Get samples of beginning and end of block
bkStart_samps = 1+find( diff(SoundData(8,:)==ib) ==  1 );
bkStop_samps  =   find( diff(SoundData(8,:)==ib) == -1 );

% Remove all blocks that are not this spl, lp noise, am depth,
% OR that contain artifact 
% OR that animal wasnt on the spout
% OR are too short in duration (when recording cut off early)
rm_bk = [];
for it = 1:numel(bkStart_samps)
    if (SoundData(8,bkStart_samps(it)-1)==11)
        rm_bk = [rm_bk it];
    end
    if any(intersect(bkStart_samps(it):bkStop_samps(it),ArtifactFlag))...
            || ~all(SoundData(3, bkStart_samps(it):bkStop_samps(it)) ==amd)...
            || ~all(SoundData(4, bkStart_samps(it):bkStop_samps(it)) ==spl)...
            || ~all(SoundData(6, bkStart_samps(it):bkStop_samps(it)) ==1)...
            || ((bkStop_samps(it)-bkStart_samps(it))/fs) < 1.9
        rm_bk = [rm_bk it];
    end
end
bkStart_samps(rm_bk) = [];
bkStop_samps(rm_bk) = [];

% Convert to ms
bkStart_ms = round( bkStart_samps / fs*1000 );
bkStop_ms  = round( bkStop_samps  / fs*1000 );

if ~isempty(bkStop_ms)
    if any((bkStop_ms-bkStart_ms)>2000)
                
        vs = find((bkStop_ms-bkStart_ms)>2000);
        
        skip_this=[];
        add_this=[];
        add_this_ms=[];
        
        for iv=1:numel(vs)
            
            % Find boundary to split in half
            newStart_samps = bkStart_samps(vs(iv)) + round((bkStop_samps(vs(iv))-bkStart_samps(vs(iv)))/2);
            newStop_samps  = bkStart_samps(vs(iv)) + round((bkStop_samps(vs(iv))-bkStart_samps(vs(iv)))/2) -1;
            
%             
%             % Plot it to see if it makes sense
%             hfv=figure;
%             set(hfv,'Position',figsize)
%             plot((bkStart_samps(vs(iv)):bkStop_samps(vs(iv)))-(bkStop_samps(vs(iv))-bkStart_samps(vs(iv))),...
%                 SoundData(2, (bkStart_samps(vs(iv)):bkStop_samps(vs(iv)))-(bkStop_samps(vs(iv))-bkStart_samps(vs(iv)))),'k')
%             hold on
%             plot(bkStop_samps(vs(iv))+(1:(bkStop_samps(vs(iv))-bkStart_samps(vs(iv)))),...
%                 SoundData(2, bkStop_samps(vs(iv))+(1:(bkStop_samps(vs(iv))-bkStart_samps(vs(iv)))) ) ,'k')            
%             
%             plot(bkStart_samps(vs(iv)):bkStop_samps(vs(iv)),SoundData(2,bkStart_samps(vs(iv)):bkStop_samps(vs(iv))))
%             
%             plot([newStart_samps newStart_samps],[-5 5])
%             xlim([bkStart_samps(vs(iv)) bkStop_samps(vs(iv))])
%             
%             
%             UserDecision = input('What should be done with this block? \nType ''s'' to skip, or ''add'' to add\n','s');
%             
%             switch UserDecision
%                 case 's'
%                     skip_this = [skip_this vs(iv)];
%                 case 'add'
%                     add_this  = [add_this newStart_samps];
%                     add_this_ms = [add_this_ms round( newStart_samps / fs*1000 )];
%             end
%             
%             close(hfv)
            
            %% auto split
            add_this  = [add_this newStart_samps];
            add_this_ms = [add_this_ms round( newStart_samps / fs*1000 )];
        end
        
        % Now make the corrections
        
        % Remove the indices of the ambiguous blocks
        bkStart_samps(skip_this) = [];
        bkStop_samps(skip_this)  = [];
        bkStart_ms(skip_this)    = [];
        bkStop_ms(skip_this)     = [];
        % Add indices of blocks to be split in half
        bkStart_samps = sort([ bkStart_samps add_this  ]);
        bkStop_samps  = sort([ bkStop_samps add_this-1 ]);
        bkStart_ms    = sort([ bkStart_ms add_this_ms  ]);
        bkStop_ms     = sort([ bkStop_ms add_this_ms-1 ]);
        
    end
end

if numel(bkStart_ms) ~= numel(bkStop_ms)
    keyboard
end


end