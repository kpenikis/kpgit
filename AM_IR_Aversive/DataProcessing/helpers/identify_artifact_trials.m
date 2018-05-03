function Info = identify_artifact_trials(Info,Phys,TrialData)
% [Info,Phys,TrialData] = identify_artifact_trials(Info,Phys,TrialData)
%  Find the trials that contain artifact. Mostly automatic, but may ask for
%  user input about borderline cases. Usually called by pp_processPhys, but
%  can also be run on partially processed sessions, by using
%  coldcall_identify_artifact. 
%  KP, 2018-02
%

window_size = 0.1; %s
window_size = round(window_size*Info.fs); %samples

ArtifactFlags_trialch = zeros(size(TrialData,1),size(Phys,1));

for ich = 1:size(Phys,1)
    if ich==8, continue, end
    
    % First, calculate rms distribution for this channel  
    fprintf('calculating RMS distribution for ch %i...',ich)
    tic 
    rms_ch = envelope(Phys(ich,:),window_size,'rms');
    fprintf('done (%0.1f s)\n',toc)
    
    rms_med = median(reshape(rms_ch,[size(rms_ch,1)*size(rms_ch,2),1]) ,'omitnan');
    rms_std = std(   reshape(rms_ch,[size(rms_ch,1)*size(rms_ch,2),1]) ,'omitnan');
    
    
    % Find where RMS is above threshold
    std_factor = 3;
    artifact_flag_rms = rms_ch > (rms_med+std_factor*rms_std);
    
    % Downsample to 1 ms
    clear artifact_flag_rms_ms
    artifact_flag_rms_ms = find(round(resample( double(artifact_flag_rms(1:10:length(artifact_flag_rms))) ,10000, round(Info.fs) )));
    
    
    %% Flag the Trials that have artifact
    
    for it = 1:size(TrialData,1)
        
        if sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))>9 && it>1
            ArtifactFlags_trialch(it,ich) = 1;
        elseif sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))<9 && sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))>0 && it>1
            % Plot borderline cases to check manually
            ArtifactFlags_trialch(it,ich) = manually_check_this_trial( it, ...
                Phys(ich, round(TrialData.onset(it)/1000*Info.fs):round(TrialData.offset(it)/1000*Info.fs)) );
        elseif it==1 && (sum(ismember(TrialData.onset(it):TrialData.offset(it),artifact_flag_rms_ms))/(TrialData.offset(it)-TrialData.onset(it)))>0.05
            % Plot trial 1 to check manually if more than 5% is flagged
            ArtifactFlags_trialch(it,ich) = manually_check_this_trial( it, ...
                Phys(ich, round(TrialData.onset(it)/1000*Info.fs):round(TrialData.offset(it)/1000*Info.fs)) );
        end
        
    end
    
end

%% Save the artifact flags into the Info struct

for ich = 1:size(Phys,1)
Info.artifact(ich).trials = find(ArtifactFlags_trialch(:,ich));
end


end




