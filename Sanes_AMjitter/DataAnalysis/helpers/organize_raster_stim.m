function stimdata = organize_raster_stim(raster)
% stimdata = organize_raster_stim(raster)
%   Separate raster entries based on the stimulus parameters. Organize by
%   block and acoustic pars, determine which stim parameter can be 
%   discriminated, then store those raster entries together. Also groups
%   raster data across blocks.
%

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


% * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * *
% Get blocks for loop and designate which ones to combine
[blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));
% * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * *


for ib = blocks
    
    % Get data for these blocks
    bk_raster = raster([raster.block]==ib);
    
    % Collapse blocks if needed
    if any(group_blocks(:,1)==ib)
        bk_raster = raster_collapse_blocks(bk_raster, raster([raster.block]== group_blocks(group_blocks(:,1)==ib,2)) );
    end
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    for ip = 1:max(np)
        
        % Get raster of this set of params, and find unique stimuli
        param_raster = bk_raster(np==ip);
        
        % Find what stimulus parameter was discriminated in this block
        JittDpthBeh_stim = unique([[param_raster.jitter]' num2str([param_raster.AMdepth]') [param_raster.behaving]'],'rows','stable');
        JittDpthBeh_stim = strtrim(JittDpthBeh_stim);
        
        discrVar=[];
        %%% xJITTER only
        if     numel(unique(JittDpthBeh_stim(:,1)))>3 ...
            && any(strcmp(JittDpthBeh_stim(:,1),'0')) ...
            && numel(unique(JittDpthBeh_stim(:,2)))<4
                discrVar = 'nm_jitter';
                
        %%% xDEPTH only
        elseif numel(unique(JittDpthBeh_stim(:,2)))>3 ...
            && any(strcmp(JittDpthBeh_stim(:,2),'0')) ...
            && numel(unique(JittDpthBeh_stim(:,1)))<4
                discrVar = 'nm_depth';

        %%% JITTER and DEPTH x axis
        elseif numel(unique(JittDpthBeh_stim(:,2)))>3 ...
            && any(strcmp(JittDpthBeh_stim(:,2),'0')) ...
            && numel(unique(JittDpthBeh_stim(:,1)))>3 ...
            && any(strcmp(JittDpthBeh_stim(:,1),'0')) 
                keyboard
                discrVar = 'nm_jitter';
                discrVar = 'nm_depth';
        
        else % e.g. increasing decreasing jitter
            keyboard
            discrVar = 'inc_dec';
            
        end
        
        
        % Make reference structure 
        stimdata(blocks==ib).block = ib;
        stimdata(blocks==ib).pars(ip).pars = LP_HP_dB_rate(ip,:);
        stimdata(blocks==ib).pars(ip).stimvals = JittDpthBeh_stim;
        stimdata(blocks==ib).pars(ip).(discrVar) = param_raster;
        
        
    end
end

end

