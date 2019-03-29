function [data,Info,TrialData] = identify_artifact_regions_KS(data,Info,TrialData)
% ABANDONED IN FAVOR OF identify_artifact_trials_64ch

ds = round(Info.fs/Info.fs_sound);


%%% Artifact identification based on RMS would work better on filtered
%%% data. Too many false positives with raw data. 

sampwin  = round(Info.fs/5);  %200 ms
sampstep = round(Info.fs/10); %100 ms
begwins  = 1:sampstep:size(dataFilt,2);
Channels = 51;

for ch = Channels
    
    %%% Artifact identification based on RMS would work better on filtered
    %%% data. Too many false positives with raw data.
    clear dataFilt
    dataFilt = filter_data( double(data(ch,:)), Info.fs );
    
    % Calculate RMS
    roughRMS = nan(size(begwins));
    for iw = 1:numel(begwins)
        roughRMS(iw) = rms(dataFilt(1,begwins(iw):min(begwins(iw)+sampwin-1,size(dataFilt,2))));
    end
    
    threshRMS = 1.3*median(roughRMS);
    
    
    % Create logical array of artifact flags
    artifact_flag_ds = zeros(size(roughRMS));
    artifact_flag_ds(roughRMS>threshRMS) = 1;
    artifact_flag = [];
    artifact_flag = logical(repmat(artifact_flag_ds,sampstep,1));
    artifact_flag = reshape(artifact_flag,1,size(artifact_flag,1)*size(artifact_flag,2));
    if size(artifact_flag,2) ~= size(dataFilt,2)
        artifact_flag = artifact_flag(1,1:size(dataFilt,2));
    end
    
    % Plot to check
    hfart=figure(1); clf
    plot(find(artifact_flag),data(ch,artifact_flag),'r','LineWidth',1.5)
    hold on
    plot(dataFilt(1,:),'k')
    keyboard
    close(hfart)
    
    
    % Find onset of drinking to denote official start of session
%     plot([SpoutTTL(1)*ds SpoutTTL(1)*ds],[-1 1].*0.1,'k','LineWidth',2)
        
    % Create white noise to match channel and replace artifact in Phys with this noise
    rms_sess_avg = rms(dataFilt(1,artifact_flag==0));
    dataFilt(1,artifact_flag) = rms_sess_avg*randn(1,sum(artifact_flag))-rms_sess_avg/2;
    
    % Also remove Phys data before actual session begins
%     data(ich,1:SpoutTTL(1)*ds-1) = 0;
    
    % View final cleaned data
    plot(dataFilt,'b')
    keyboard
    %if good, proceed to replace data with cleaned data
    data(ch,:) = dataFilt;
    
    
    % Downsample the artifact flags to match sampling rate of TrialData
    artifact_flag = artifact_flag(1,1:ds:size(data,2));
    
    if size(artifact_data,2) ~= size(TrialData,2)
        artifact_data = artifact_data(1,1:size(TrialData,2));
    end
    
    % Save the artifact data for this channel
    Info.artifact(ich).SDsamples  = find(artifact_data);
    
    
end %ich


%% Then launch gui to flag entire stimulus blocks

% [Info,TrialData] = manually_check_artifact_blocks(Info,data,TrialData,artifact_samples);





end




