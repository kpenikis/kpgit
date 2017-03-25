function [stimdata,clu_raster] = get_raster_indices(stimdata,clu_raster,thisblock,indVar)

for id = 1:numel(stimdata.stimvals)
    
    thisdepth = stimdata.stimvals(id).depth;
    
    raster_idxs=[];
    for ij = 1:numel(stimdata.stimvals(id).jitter)
        
        thisjitter = stimdata.stimvals(id).jitter(ij);
        
        % Find raster struct entries that match this depth and jitter
        raster_idx = intersect(find(strcmp([clu_raster.jitter],thisjitter)),find([clu_raster.AMdepth]==thisdepth));
        
        % Find matching behavioral state
        raster_idx(~strcmp({clu_raster(raster_idx).behaving},stimdata.stimvals(id).behav{ij})) = [];
        
        % Remove entries whos pars dont match
        raster_idx([clu_raster(raster_idx).HP] ~= stimdata.pars(1)) = [];
        raster_idx([clu_raster(raster_idx).LP] ~= stimdata.pars(2)) = [];
        raster_idx([clu_raster(raster_idx).dB] ~= stimdata.pars(3)) = [];
        raster_idx([clu_raster(raster_idx).AMrate] ~= stimdata.pars(4)) = [];
        
        % if data from two blocks are supposed to be combined
        %%%%% CHECK IF THERE IS AN ENTRY WITH >1 BLOCK (FOR RERUNS)
        if numel(thisblock)>1
            
            rmv_idx=raster_idx;
            for ib = 1:numel(thisblock)
                rmv_idx([clu_raster(raster_idx).block]==thisblock(ib)) = [];
            end
            if ~isempty(rmv_idx)
                raster_idx(raster_idx==rmv_idx) = [];
            end
            
            if id==1 && ij==1
                warning('adding fields to raster struct. dat cool?')
            end
            clu_raster(end+1) = collapse_blocks(clu_raster(raster_idx));
            raster_idx = numel(clu_raster);
            
        else %typical case
            
            raster_idx([clu_raster(raster_idx).block]~=thisblock) = [];
            
        end
        
        if ~(numel(raster_idx)==1)
            % There SHOULD be one index remaining, if identical blocks
            % were combined during initial analysis phase, and nothing
            % else went wrong.
            keyboard
        end
        
        raster_idxs = [raster_idxs raster_idx];
        
    end
    
    % Add field to stim struct of cludata
    stimdata.stimvals(id).raster_idxs = raster_idxs;
    
end

end