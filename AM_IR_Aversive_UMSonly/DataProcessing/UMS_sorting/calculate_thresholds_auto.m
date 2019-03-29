function [thresh, reject] = calculate_thresholds_auto( Phys, TrialData, Info )
%  [thresh, reject] = calculate_thresholds( Phys, trial_length_s )
%    Called by pp_sort_session. Launches a gui in which a random trial is
%    selected and filtered Phys data is displayed for 4 channels at a time.
%    Waits for user input, to designate trial as clean or to skip to next
%    random trial. Once a sufficient amount of data has been designated
%    "clean" for each channel, thresholds for spike event detection and
%    artifact rejection are calculated based on the std of data, for each
%    channel.
%  
%  KP, 2016-04; last updated 2018-02
%  2018-02: now written for TrialData instead of SoundData format
% 



try
% Set some initial params
totalms = 16000;
msgathered  = 0;

% Create empty vectors
std_cln = zeros(1,size(Phys,1));
clean_samples = [];

% Cycle through randomly chosen trials
rng('shuffle')
for it = 1+randperm(size(TrialData,1)-2) %skipping trial #1 (silence) and last
    
    % Check this trial for artifact
    use_trial = 1; %default is to use it, mark with 0 if noisy
    for ich = 1:size(Phys,1)
        if ~isempty( intersect( Info.artifact(ich).trials, it ) )
            use_trial = 0;
        end
    end
    
    % Only if all channels are clean, save this trial
    if use_trial == 1
        
        clean_samples = [clean_samples round(TrialData.onset(it)/1000*Info.fs):round(TrialData.offset(it)/1000*Info.fs)];
        
        msgathered = msgathered + round(TrialData.offset(it)-TrialData.onset(it));
        if msgathered < totalms
            continue
        else
            fprintf('  %i ms of clean data gathered\n',msgathered)
            break
        end
    end
    
end %it


%% Calculate std and thresholds for each channel
for ich = 1:size(Phys,1)
    std_cln(1,ich) = mean(std( Phys(ich,clean_samples) ,0,2));
end


catch
    keyboard
end

%% Calculate final threshold values

thresh = std_cln .*  -4 ;
reject = std_cln .* -20 ;


%% Check that thresholds make sense by plotting with clean traces
% for ich = 1
%     f=figure;
%     scrsz = get(0,'ScreenSize');
%     set(f,'Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
%     for isp = 1:2
%         plot( Phys(ich-1+isp,:) + (isp*8)*10^-4 ,'Color',[0.4 0.4 0.4])
%         hold on
%         plot([0 size(Phys,2)], [(thresh(ich-1+isp) + (isp*8)*10^-4) (thresh(ich-1+isp) + (isp*8)*10^-4)], ':', 'Color', [0 0.7 0], 'LineWidth',3)
%         plot([0 size(Phys,2)], [(reject(ich-1+isp) + (isp*8)*10^-4) (reject(ich-1+isp) + (isp*8)*10^-4)], ':', 'Color', [0.7 0 0], 'LineWidth',3)
% 
%     end    
%     title([ 'Check thresholds for channels ' num2str(ich) ' & ' num2str(ich+1) ])
% end

end
