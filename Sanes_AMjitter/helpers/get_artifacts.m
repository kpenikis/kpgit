function artifact = get_artifacts(data_raw, fs)
% called by pp_prepare_format


window_size = 0.5; %s
window_size = round(window_size*fs); %samples

% af = figure;

for ich = 1:size(data_raw,3)
    if ich==8, continue, end
    
    % Calculate rms for this channel  
    data_ch = data_raw(:,:,ich);
    rms_tr = nan(size(data_raw,1), ceil(size(data_raw,2)/window_size));
    
    for itr = 1:size(data_raw,1)
        nseg=0;
        for iseg = 1:window_size:round(size(data_raw,2))
            nseg=nseg+1;
            rms_tr(itr,nseg) = sqrt(mean(data_raw(itr, iseg:(min(iseg+window_size-1,end)) ,ich).^2));
        end
    end
    rms_ch  = median(reshape(rms_tr,[size(rms_tr,1)*size(rms_tr,2),1]),'omitnan');
    rms_std = std(reshape(rms_tr,[size(rms_tr,1)*size(rms_tr,2),1]),'omitnan');
        
    amp_std(ich) = median(std( data_raw( :,:,ich ) ,0,2) ,'omitnan') ;

    
    a_trs=[];
    for itr = 1:size(data_raw,1)
        
        data = data_raw(itr,:,ich);
        
        % RMS above threshold?
        std_factor = 10;
        artifact_flag_amp = any( data>(amp_std(ich)*std_factor) | data<(-amp_std(ich)*std_factor));
        artifact_flag_rms = any(rms_tr(itr,:) > (rms_ch+std_factor*rms_std));
        
        if artifact_flag_rms
            a_trs=[a_trs itr];
        end
        
% % %         % Plot
% % %         figure(af); clf
% % %         if artifact_flag_rms
% % %             plot(data,'r')
% % %         else
% % %             plot(data,'k')
% % %         end
% % %         hold on
% % %         plot( [0 size(data_raw,2)], [(amp_std(ich)*std_factor) (amp_std(ich)*std_factor)], 'r')
% % %         plot( [0 size(data_raw,2)], [-(amp_std(ich)*std_factor) -(amp_std(ich)*std_factor)], 'r')
% % % %         xlim([0 size(data_raw,2)])
% % % %         ylim([-amp_std(ich)*30 amp_std(ich)*25])
% % %         title(['ch ' num2str(ich) ', tr ' num2str(itr)])
% % % %         
% % %         pause(1)
                        
    end
    artifact(ich).trials = a_trs;
    
end



end