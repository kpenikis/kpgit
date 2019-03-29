
function FLAG_OUT = manually_check_this_trial(it,data)
% called by identify_artifact_trials
% 
% KP, 2018-02
%


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


FLAG_OUT = 0;
user_choice = 0;

% Plot data from this flagged trial
plot( data, 'Color',[0.7 0.7 0.7])

set(gca, 'xlim', [1 length(data)],...
    'ylim',[1.5*min(data) 1.6*max(data)],...
    'XTick',[], 'YTick',[])
title([ ' Is trial ' num2str(it) ' clean or dirty?'])

% Wait for user to mark this flagged trial as clean or dirty
while user_choice==0
    pause(0.5);
end

switch user_choice
    case 1 %clean
        FLAG_OUT = 0;
    case 2 % dirty
        FLAG_OUT = 1;
end

close(hf)

end
