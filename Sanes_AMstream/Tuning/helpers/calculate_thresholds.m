function [thresh, reject] = calculate_thresholds( Phys, Info, Stim )
%  [thresh, reject] = calculate_thresholds( Phys, trial_length_s )
%    Called by pp_sort_session. Launches a gui in which a random trial is
%    selected and filtered Phys data is displayed for 4 channels at a time.
%    Waits for user input, to designate trial as clean or to skip to next
%    random trial. Once a sufficient amount of data has been designated
%    "clean" for each channel, thresholds for spike event detection and
%    artifact rejection are calculated based on the std of data, for each
%    channel.
%  
%   >> Now runs automatically, leaving out flagged trials.
%  
%  KP, 2016-04; last updated 2017-06
% 


all_clean=[];
ntrials = floor(16 / (diff(Info.t_win_ms)/1000));   % need ~16 s of clean data

for channel = 1:size(Phys,3)
    
    clean_trials = nan(1,ntrials);
    nchosen = 0;
    
    flagged = Info.artifact_trs(channel).manual;
    
    
    % Cycle through trials, skipping flagged ones
    
    rng('shuffle')
    for it = randperm(size(Stim,2)-1)
        
        if isempty(intersect(it,flagged))
            nchosen = nchosen+1;
            clean_trials(:,nchosen) = it;
        end
        
        if nchosen==ntrials
            break
        end
        
    end
    
    all_clean = [all_clean; clean_trials];
    
end


% Calculate std and thresholds for each trial
for ich = 1:size(Phys,3)
    std_cln(1,ich) = mean(std(Phys( all_clean(ich,:) ,:,ich),0,2));
        
end

thresh = std_cln .*  -4 ;
reject = std_cln .* -20 ;


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









