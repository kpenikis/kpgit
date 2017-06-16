function [Info,SoundData] = identify_artifact_regions(Info,Phys,SoundData)


ds = round(Info.fs/Info.fs_sound);

window_size = 0.1; %s
window_size = round(window_size*Info.fs); %samples

artifact_samples = nan(size(Phys,1),size(SoundData,2));

% af = figure;

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
    
    
    %%%% DOWNSAMPLE TO 1 MS
    
    artifact_data = repmat(artifact_flag_rms,window_size,1);
    artifact_data = reshape(artifact_data,1,size(artifact_data,1)*size(artifact_data,2));
    artifact_data = artifact_data(1,1:ds:size(Phys,2));
    
    if size(artifact_data,2) ~= size(SoundData,2)
        artifact_data = artifact_data(1,1:size(SoundData,2));
    end
    
    % Save this data for each channel
    artifact_samples(ich,:) = artifact_data;
    
    
    % Plot
% %     figure; clf
% %     plot(Phys(ich,1:ds:end),'k')
% %     hold on
% %     plot(artifact_data.*(10^-3),'r')
% %     
% %     figure; clf
% %     plot(rms_ch,'k')
% %     hold on
% %     plot(artifact_flag_rms.*10^-4,'r')
% %     plot([1 length(rms_ch)],[(rms_med+std_factor*rms_std) (rms_med+std_factor*rms_std)],'b')
% %             
% %     xlim([1 size(rms_ch,2)])
% %     title(['ch ' num2str(ich)])
% %     pause(3)
    
end


%% Then launch gui to flag entire stimulus blocks

[Info,SoundData] = manually_check_artifact_blocks(Info,Phys,SoundData,artifact_samples);





end




