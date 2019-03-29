function [Info,Phys,SoundData] = identify_artifact_regions_KS(Info,Phys,SoundData)


ds = round(Info.fs/Info.fs_sound);

window_size = 0.1; %s
window_size = round(window_size*Info.fs); %samples


for ich = 1:size(Phys,1)
    if ich==8, continue, end
    
    % First, calculate rms distribution for this channel  
    rms_ch = nan( 1, ceil(size(Phys,2)/window_size) );
    
    nseg=0;
    for iseg = 1:window_size:size(Phys,2)
        nseg=nseg+1;
        rms_ch(1,nseg) = sqrt(mean( Phys(ich, iseg:(min(iseg+window_size-1,end)) ).^2 ));
    end
    
    rms_med = median(reshape(rms_ch,[size(rms_ch,1)*size(rms_ch,2),1]) ,'omitnan');
    rms_std = std(   reshape(rms_ch,[size(rms_ch,1)*size(rms_ch,2),1]) ,'omitnan');
        
    amp_std(ich) = median(std( Phys( ich,: ) ,0,2) ,'omitnan') ;

    
    % RMS above threshold?
    std_factor = 3;
    artifact_flag_amp = rms_ch > (amp_std(ich)*std_factor) | rms_ch < (-amp_std(ich)*std_factor);
    artifact_flag_rms = rms_ch > (rms_med+std_factor*rms_std);
    
    artifact_data = [];
    artifact_data = repmat(artifact_flag_rms,window_size,1);
    artifact_data = reshape(artifact_data,1,size(artifact_data,1)*size(artifact_data,2));
    
    if size(artifact_data,2) ~= size(Phys,2)
        artifact_data = artifact_data(1,1:size(Phys,2));
    end
    
    
    % Plot
%     figure(1); clf
%     plot(Phys(ich,:),'k')
%     hold on
%     plot(find(artifact_data),Phys(ich,artifact_data),'r')
    
    % Find onset of drinking to denote official start of session
    SpoutTTL = find(SoundData(7,:));
%     plot([SpoutTTL(1)*ds SpoutTTL(1)*ds],[-1 1].*0.005,'k','LineWidth',2)
        
    % Create white noise to match channel and replace artifact in Phys with this noise
    rms_sess_avg = rms(Phys(ich,~artifact_data));
    Phys(ich,artifact_data) = rms_sess_avg*randn(1,sum(artifact_data))-rms_sess_avg/2;
    
    % Also remove Phys data before actual session begins
    Phys(ich,1:SpoutTTL(1)*ds-1) = 0;
    
    % View final cleaned data
%     plot(Phys(ich,:),'b')
    
    
    % Downsample the artifact flags to match sampling rate of SoundData
    artifact_data = artifact_data(1,1:ds:size(Phys,2));
    
    if size(artifact_data,2) ~= size(SoundData,2)
        artifact_data = artifact_data(1,1:size(SoundData,2));
    end
    
    % Save the artifact data for this channel
    Info.artifact(ich).SDsamples  = find(artifact_data);
    
    
end %ich


%% Then launch gui to flag entire stimulus blocks

% [Info,SoundData] = manually_check_artifact_blocks(Info,Phys,SoundData,artifact_samples);





end




