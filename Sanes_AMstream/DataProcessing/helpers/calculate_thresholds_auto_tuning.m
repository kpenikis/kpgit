function [thresh, reject] = calculate_thresholds_auto_tuning( Phys, SoundData, Info )
%  [thresh, reject] = calculate_thresholds( Phys, trial_length_s )
%    Called by pp_sort_session. Launches a gui in which a random trial is
%    selected and filtered Phys data is displayed for 4 channels at a time.
%    Waits for user input, to designate trial as clean or to skip to next
%    random trial. Once a sufficient amount of data has been designated
%    "clean" for each channel, thresholds for spike event detection and
%    artifact rejection are calculated based on the std of data, for each
%    channel.
%  
%  KP, 2016-04; last updated 2017-06
% 


% Set some initial params
nblocks = floor(16 / 2);   % need ~16 s of clean data
ds = round(Info.fs/Info.fs_sound);

% Find the indices of the start of blocks
AllBlocks = 1+find(diff([SoundData(8,:) 0]));

% Create empty vector for each channel std
std_cln = zeros(1,size(Phys,1));

for channel = [1 5 9 13]
    
    clean_samples = [];
    nchosen = 0;
    
    % Cycle through randomly chosen segments, user selects clean ones
    rng('shuffle')
    for ibs = randperm(size(AllBlocks,2))
        
        % Skip if unmodulated noise, silence, or irrelevant time in
        % recording
        try
        if SoundData(8,AllBlocks(ibs)) > 10 || SoundData(8,AllBlocks(ibs))==0
            continue
        end
        catch
            keyboard
        end
        
        % Check these segments for artifact
        use_segment = 1; %default is to use it, mark with 0 if noisy
        for isp = 0:3
            if ~isempty( intersect( Info.artifact(channel+isp).SDsamples, AllBlocks(ibs):(AllBlocks(ibs+1)-1) ) )
                use_segment = 0;
            end
        end
        
        % Only if all 4 channels are clean, save these samples
        if use_segment == 1
            
            clean_samples = [clean_samples AllBlocks(ibs):(AllBlocks(ibs+1)-1)];
            
            nchosen = nchosen+1;
            if nchosen<nblocks
                continue
            else
                fprintf('  %i clean segments gathered for chs %i-%i\n',nblocks, channel, channel+3)
                break
            end
        end
        
    end
    
    % Calculate std and thresholds for each channel
    for isp = 0:3
        std_cln(1,channel+isp) = mean(std( Phys(channel+isp,clean_samples) ,0,2));
    end
        
end

% Calculate final threshold values
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




