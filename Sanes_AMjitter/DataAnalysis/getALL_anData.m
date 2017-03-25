function getALL_anData(subject,session)
% getALL_anData(subject,session)
%   Plot responses of all clusters of a session to all stimuli, for all
%   response metrics.
%
%   Inputs
%     subject: subject name as string
%     session: session label as string
%     channel: (optional) channel number as double
%         clu: (optional) cluster label as double
%
% KP, 02-2017,
%

global fn

% Load Data structure
fn = set_paths_directories(subject,session);
filename = sprintf('%s_sess-%s_Data',subject,session);
load(fullfile(fn.processed,subject,filename));

% Go through each cluster and call analysis programs

allclusters = fieldnames(Data);
for unit = allclusters'
    
    channel = Data.(unit{:}).labels(1,3);
    clu     = Data.(unit{:}).labels(1,1);
    
    run_analyses(Data.(unit{:}).raster, subject,session,channel,clu,unit{:})
    
    close all
    
end


end


function run_analyses(raster, subject,session,channel,clu,cluname)

global fn

METRICS = {'FR' 'FF' 'FF-avPds' 'FF-Pds' 'VS' 'VS-Pds' 'RS' 'RS-Pds' 'standardFR' 'Corr'};

cludata = struct;

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

% Get blocks for loop and designate which ones to combine
[blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));

for ib = blocks
    
    % Get data for these blocks
    bk_raster = raster([raster.block]==ib);
    
    % Add data from blocks set to combine, if needed
    if ~isempty(group_blocks(:,1)==ib)
        bk_raster = [bk_raster raster([raster.block]== group_blocks(group_blocks(:,1)==ib,2) )];
        bk_str = [num2str(ib) num2str(group_blocks(group_blocks(:,1)==ib,2))];
        bk_vec = [ib group_blocks(group_blocks(:,1)==ib,2)];
    else
        bk_str = num2str(ib);
        bk_vec = ib;
    end
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    for ip = 1:max(np)
        
        stimvals=struct;
        
        % Make data structure for this cluster
        cludata.blk(blocks==ib).block = bk_vec;
        cludata.blk(blocks==ib).pars(ip).pars = LP_HP_dB_rate(ip,:);
        cludata.blk(blocks==ib).pars(ip).stimvals = struct;
        
        % Get raster of this set of params, and find unique stimuli
        param_raster = bk_raster(np==ip);
        %                 [depths,~,ndpth] = unique([param_raster.AMdepth]);
        %                 [jitters,~,njit] = unique([param_raster.jitter],'stable');
        %                 %%% xJITTER only
        %                 if( numel(depths(depths~=0))<2 && numel(jitters)>1 )
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for im = 1:numel(METRICS)
            
            [datastruct,fieldname,stimvals] = an_xJitter(subject,session,channel,clu,...
                METRICS{im},param_raster,stimvals);
            
            cludata.blk(blocks==ib).pars(ip).(fieldname) = datastruct;
        end
        cludata.blk(blocks==ib).pars(ip).stimvals = stimvals;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %                     %%% xDEPTH only
        %                 elseif( numel(depths(depths~=0))>1 && numel(jitters)<2 )
        %                     %%% JITTER and DEPTH x axis
        %                 elseif( numel(depths(depths~=0))>1 && numel(jitters)>1 )
        %                 end
    end
end


% Save to cluster data file
filename = sprintf( '%s_sess-%s_%s',subject,session,cluname);
load(fullfile(fn.sess_data,filename))
eval(sprintf('%s.block = cludata.blk;',cluname));
save(fullfile(fn.sess_data,filename),cluname,'-v7.3');



end % function runALL_anPlots


