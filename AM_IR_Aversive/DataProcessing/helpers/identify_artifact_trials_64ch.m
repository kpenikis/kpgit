function [Info,DATA] = identify_artifact_trials_64ch(Info,DATA,TrialData)
% [Info,DATA] = identify_artifact_trials(Info,DATA,TrialData)
%  Find the trials that contain artifact. Mostly automatic, but may ask for
%  user input about borderline cases. Usually called by pp_processPhys, but
%  can also be run on partially processed sessions, by using
%  coldcall_identify_artifact. 
%  KP, 2019-01
%


Shanks = unique(Info.probeMap.kcoords(Info.probeMap.connected));

window_size = 0.1; %s
window_size = round(window_size*Info.fs); %samples

ArtifactFlags_TrGrp = zeros(size(TrialData,1),max(Shanks));

% To save time, calculate artifact for one whole channel group (/shank) at a time
for ishank = Shanks
    
    % Select random chnnel from this group
    channels = Info.probeMap.chanMap( Info.probeMap.kcoords==ishank & Info.probeMap.connected'==1);
    foo = randperm(length(channels));
    ich = channels(foo(1));
    
    
    % First, filter the channel
    %%% (artifact identification based on RMS yields many false positives
    %%% using unfiltered raw data)
    clear dataFilt
    fprintf('  filtering group %i, ch %i... ',ishank,ich)
    tic 
    dataFilt = filter_data( double(DATA(ich,:)), Info.fs );
    fprintf('done (%0.1f s)\n',toc)
    
    
    % Calculate rms distribution for this channel  
    fprintf('  calculating RMS distribution for group %i, ch %i... ',ishank,ich)
    tic 
    rms_ch = envelope(single(dataFilt),window_size,'rms');
    fprintf('done (%0.1f s)\n',toc)
    
    
    % Find where RMS is above threshold
    std_factor  = 1;
    thresh1 = median(rms_ch) + std(rms_ch) * std_factor;
    mult_factor = 1.5;
    thresh2 = median(rms_ch) * mult_factor;
    
    artifact_flag_rms = rms_ch > thresh2;
    
    
    % Check while debugging
%     keyboard
    
%     hfart = figure(1); clf
%     plot(find(artifact_flag_rms),DATA(ich,artifact_flag_rms),'r','LineWidth',1.5)
%     hold on
%     plot(dataFilt,'k')
%     
%     hfart = figure(1); clf
%     plot(rms_ch,'k'); hold on
%     plot([0 length(dataFilt)],[thresh2 thresh2],'b')
%     
%     hfart = figure(1); clf
%     plot(find(artifact_flag_rms),DATA(ich,artifact_flag_rms),'r','LineWidth',1.5)
%     hold on
%     plot(DATA(54,:),'k')
%     
%     hfart = figure(1); clf
%     plot(find(artifact_flag_rms),rms_ch(1,artifact_flag_rms),'r.','LineWidth',1.5)
%     hold on
%     for iich = channels(2:2:8)
%         clear dataFilt rms_ch
%         dataFilt = filter_data( double(DATA(iich,:)), Info.fs );
%         rms_ch = envelope(single(dataFilt),window_size,'rms');
%         plot(rms_ch)
%     end
%     
%     close(hfart)
    
    
    % Downsample to 1 ms
    clear artifact_flag_rms_ms
    artifact_flag_rms_ms = find(round(resample( double(artifact_flag_rms(1:10:length(artifact_flag_rms))) ,10000, round(Info.fs) )));
    
    
    %% Flag the Trials that have artifact
% %     
% %     ms_thresh = 50; % 5 percent of pdc trial duration
% %     
% %     for it = 1:size(TrialData,1)
% %         
% %         if sum( ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms) )>ms_thresh && it>1
% %             % More than 50 ms of artifact, flag it
% %             ArtifactFlags_TrGrp(it,ishank) = 1;
% %             
% %         elseif sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))<ms_thresh && sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))>0 && it>1
% %             % Plot borderline cases (less than 50 ms) to check manually 
% %             ArtifactFlags_TrGrp(it,ishank) = manually_check_this_trial( it, ...
% %                 dataFilt(1, round(TrialData.onset(it)/1000*Info.fs):round(TrialData.offset(it)/1000*Info.fs)) );
% %         
% %         elseif it==1 && (sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))/(TrialData.offset(it)-TrialData.onset(it)))>0.05
% %             % Plot trial 1 to check manually if more than 5 percent is flagged
% %             ArtifactFlags_TrGrp(it,ishank) = manually_check_this_trial( it, ...
% %                 dataFilt(1, round(TrialData.onset(it)/1000*Info.fs):round(TrialData.offset(it)/1000*Info.fs)) );
% %         end
% %         
% %     end %it
% %     
    
    %% Replace artifact periods with white noise, necessary to clean up sorting
    
    rms_sess_avg = rms(dataFilt(1,artifact_flag_rms==0));
    DATA(channels,artifact_flag_rms) = rms_sess_avg*randn(length(channels),sum(artifact_flag_rms))-rms_sess_avg/2;
    
    % Replace period before session starts with white noise
    DATA(channels,1:TrialData.onset(1)) = rms_sess_avg*randn(length(channels),TrialData.onset(1))-rms_sess_avg/2;
    
    
    %% Save the artifact flags to the Info struct
    
%     for ich = channels
%         Info.artifact(ich).trials = find(ArtifactFlags_TrGrp(:,ishank));
%     end
    
    
end %ishank


% % % Apply result from sample channel of group to all channels in group
% % for ich = 1:size(Phys,1)
% %     Info.artifact(ich).trials = find(ArtifactFlags_TrGrp(:,ich));
% % end


end




