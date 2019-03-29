
function Info = manually_check_artifact_blocks(Info,Phys,SoundData,artifact_samples)
% Go through stimulus blocks that have include segments flagged based on 
% RMS deviations, and mark clean or dirty.



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Prepare gui figure
set(0,'DefaultAxesFontSize',14)

% Button callback function
    function buttonCallback(newVal)
        user_choice = newVal;
    end

% Set up figure
hf=figure;
scrsz = get(0,'ScreenSize');
set(hf,'Position',[1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)/2]);
ax = axes('Parent',hf,'position',[0.05 0.2 0.9 0.75]);

% Create push buttons
cln = uicontrol('Style', 'pushbutton', 'String', 'Clean',...
    'Units', 'normalized', ...
    'Position', [0.25 0.05 0.15 0.1],...
    'BackgroundColor', [0 0.8 0], ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'Callback', @(src,evnt)buttonCallback(1) ) ;
nxt = uicontrol('Style', 'pushbutton', 'String', 'Dirty',...
    'Units', 'normalized', ...
    'Position', [0.6 0.05 0.15 0.1],...
    'BackgroundColor', [0.8 0 0], ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'Callback', @(src,evnt)buttonCallback(2)) ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ds = round(Info.fs/Info.fs_sound);
ind_endArt = [];
bbb = [];
ccc = [];
ddd = [];

AllBlocks = 1+find(diff([SoundData(9,:) 0]));
% SoundData(8, (trF(1)-1):(trF(1)) )

% Go through each channel (eventually make a more efficient way of doing
% this - 4 at a time)
for ich = 1:size(Phys,1)
        
    flaggedSamples = artifact_samples(ich,:);
    ind_startArt = find(diff(flaggedSamples)==1)+1;
    ind_endArt   = find(diff(flaggedSamples)==-1);
    markedSamples  = zeros(size(flaggedSamples));
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    user_choice = 0;
    for ifs = 1:numel(ind_startArt)

        samples = ind_startArt(ifs) : ind_endArt(ifs);
        
        % Plot data from this flagged block
        try
            plot( Phys(ich, ds.* [(samples(1)-round(Info.fs)):(samples(end)+round(Info.fs))] ) , 'Color',[0.7 0.7 0.7])
            hold on
            plot( round(Info.fs) + (1:length(samples)), Phys(ich, ds.* samples) , 'Color','r')
        catch
            keyboard
%             if size(Phys,2) < (ds*samples(end))
%                 samples = samples(1:end-1);
%             end
%             plot( round(Info.fs_sound+1):(round(Info.fs_sound)+length(samples)), Phys(ich,ds.*samples), 'Color',[0.3 0.3 0.3])
%             hold on
        end
        set(gca, 'xlim', [0 length(samples)+2*round(Info.fs)], 'XTick',[], 'YTick',[])
        title([ 'ch ' num2str(ich) ])
        hold off
        
        % Wait for user to mark this flagged trial as clean or dirty
        while user_choice==0
            pause(0.5);
        end
        
        switch user_choice
            case 1 %clean
                user_choice = 0;
                continue
            case 2 % dirty
                markedSamples(samples) = 1;
                user_choice = 0;
                continue
        end
        
    end %ifs
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Save manually marked to Info struct
%     Info.artifact(ich).flagged = find(flaggedSamples);
    Info.artifact(ich).SDsamples  = find(markedSamples);
    
    
end %ich

close(hf)

end
