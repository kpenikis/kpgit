function [thresh, reject] = calculate_thresholds( Phys, SoundData, Info )
%  [thresh, reject] = calculate_thresholds( Phys, trial_length_s )
%    Called by pp_sort_session. Launches a gui in which a random trial is
%    selected and filtered Phys data is displayed for 4 channels at a time.
%    Waits for user input, to designate trial as clean or to skip to next
%    random trial. Once a sufficient amount of data has been designated
%    "clean" for each channel, thresholds for spike event detection and
%    artifact rejection are calculated based on the std of data, for each
%    channel.
%  
%  KP, 2016-04; last updated 2017-06
% 


set(0,'DefaultAxesFontSize',14)

% Button callback function
    function buttonCallback(newVal)
        user_choice = newVal;
    end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set up figure

f=figure;
scrsz = get(0,'ScreenSize');
set(f,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
ax = axes('Parent',f,'position',[0.05 0.2 0.9 0.75]);

% Create push buttons
cln = uicontrol('Style', 'pushbutton', 'String', 'Clean',...
    'Units', 'normalized', ...
    'Position', [0.25 0.05 0.15 0.1],...
    'BackgroundColor', [0 0.8 0], ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'Callback', @(src,evnt)buttonCallback(1) ) ;

nxt = uicontrol('Style', 'pushbutton', 'String', 'Next',...
    'Units', 'normalized', ...
    'Position', [0.6 0.05 0.15 0.1],...
    'BackgroundColor', [0.8 0 0], ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'Callback', @(src,evnt)buttonCallback(2)) ;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


nblocks = floor(16 / 2);   % need ~16 s of clean data
ds = round(Info.fs/Info.fs_sound);

AllBlocks = 1+find(diff([SoundData(8,:) 0]));

std_cln = zeros(1,size(Phys,1));

for channel = [1 5 9 13]
    
    clean_samples = [];
    nchosen = 0;
    user_choice = 0;    
    
    % Cycle through randomly chosen segments, user selects clean ones
    
    rng('shuffle')
    for ibs = randperm(size(AllBlocks,2))
        
        % Skip if unmodulated noise, silence, or irrelevant time in
        % recording
        if SoundData(8,AllBlocks(ibs)) > 10 || SoundData(8,AllBlocks(ibs))==0
            continue
        end
        
        % Plot data from randomly selected trial
        cla(ax); hold on
        for isp = 1:4
            plot( Phys(channel-1+isp, ds.*(AllBlocks(ibs):(AllBlocks(ibs+1)-1)) ) + (isp*8)*10^-4 ,'Color',[0.4 0.4 0.4])
        end
        set(gca, 'XTick',[], 'YTick',[])
        title([ num2str(nchosen) ' / ' num2str(nblocks) ' clean trials chosen so far  (channels ' num2str(channel) ' - ' num2str(channel+3) ')' ])
        
        
        
        % Wait for user to skip to next channel or label it clean
        
        while user_choice==0
            pause(0.5);
        end
        
        switch user_choice
            case 1 %clean
                nchosen = nchosen+1;
%                 clean_blocks(:,nchosen) = ibs;
                clean_samples = [clean_samples AllBlocks(ibs):(AllBlocks(ibs+1)-1)];
                user_choice = 0;
                if nchosen<nblocks
                    continue
                else
                    fprintf('  %i clean segments gathered for chs %i-%i\n',nblocks, channel, channel+3)
                    break
                end
                
            case 2 % next
                user_choice = 0;
                continue
                
        end
        
    end
    
    % Calculate std and thresholds for each channel
    for isp = 0:3
        std_cln(1,channel+isp) = mean(std( Phys(channel+isp,clean_samples) ,0,2));
    end
        
end

% Calculate std and thresholds for each trial
% for ich = 1:size(Phys,3)
%     std_cln(1,ich) = mean(std(Phys( all_clean(ich,:) ,:,ich),0,2));
%         
% end

thresh = std_cln .*  -4 ;
reject = std_cln .* -20 ;

close(f);



% Check that thresholds make sense by plotting with clean traces

% for ich = 1:2:15
%     figure;
%     scrsz = get(0,'ScreenSize');
%     set(f,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
%     for isp = 1:2
%         plot_tr = all_clean(ich,randi(size(all_clean,2)));
%         plot( Phys(plot_tr,:,ich-1+isp) + (isp*8)*10^-4 ,'Color',[0.4 0.4 0.4])
%         hold on
%         plot([0 size(Phys,2)], [(thresh(ich-1+isp) + (isp*8)*10^-4) (thresh(ich-1+isp) + (isp*8)*10^-4)], ':', 'Color', [0 0.7 0], 'LineWidth',3)
%         plot([0 size(Phys,2)], [(reject(ich-1+isp) + (isp*8)*10^-4) (reject(ich-1+isp) + (isp*8)*10^-4)], ':', 'Color', [0.7 0 0], 'LineWidth',3)
% 
%     end    
%     title([ 'Check thresholds for channels ' num2str(ich) ' & ' num2str(ich+1) ])
% end
% 
% keyboard




end









