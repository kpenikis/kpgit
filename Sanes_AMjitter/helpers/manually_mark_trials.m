
function Info = manually_mark_trials(Info_old,session,Stim)
% Go through trials that have been flagged based on RMS deviations, and
% mark clear or dirty.

subject = Info_old.subject;

Info = Info_old;

flagged = Info_old.artifact_trs;

% Load Phys struct
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
fprintf('\nloading Phys struct...\n')
filename = sprintf('%s_sess-%s_Phys',subject,session);
load(fullfile(datadir,subject,filename));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Prepare gui figure
set(0,'DefaultAxesFontSize',14)

% Button callback function
    function buttonCallback(newVal)
        user_choice = newVal;
    end

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
nxt = uicontrol('Style', 'pushbutton', 'String', 'Dirty',...
    'Units', 'normalized', ...
    'Position', [0.6 0.05 0.15 0.1],...
    'BackgroundColor', [0.8 0 0], ...
    'FontSize', 24, ...
    'FontWeight', 'bold', ...
    'Callback', @(src,evnt)buttonCallback(2)) ;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Go through each channel (eventually make a more efficient way of doing
% this - 4 at a time)
for ich = 1:size(Phys,3)
    
    trF = flagged(ich).trials;
    
    trB = [];
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    user_choice = 0;
    for itr = trF
        
        % Plot data from this flagged trial
        
        %     for isp = 1:4
        %         plot( Phys(it,:,channel-1+isp) + (isp*8)*10^-4 ,'Color',[0.4 0.4 0.4])
        %     end
        plot( Phys(itr,:,ich), 'Color',[0.4 0.4 0.4])
        hold on
        set(gca, 'XTick',[], 'YTick',[])
        title([ 'ch ' num2str(ich) ' trial ' num2str(itr)])
        
        yax = get(ax,'YLim');
        plot(round([-Info.t_win_ms(1)/1000*Info.fs -Info.t_win_ms(1)/1000*Info.fs]), yax, 'b:', 'LineWidth', 2);
        plot(round([(-Info.t_win_ms(1)+Stim(itr).stimDur)/1000*Info.fs (-Info.t_win_ms(1)+Stim(itr).stimDur)/1000*Info.fs]), yax, 'b:', 'LineWidth', 2);
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
                trB = [trB itr];
                user_choice = 0;
                continue
        end
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Save manually marked to Info struct
    Info.artifact_trs(ich).manual = trB;
    
end


% Save new Info struct
filename = sprintf('%s_sess-%s_Info',subject,session);
save(fullfile(datadir,subject,filename),'Info','-v7.3');



end
