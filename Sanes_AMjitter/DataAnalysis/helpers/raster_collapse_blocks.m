function collapsed_raster = raster_collapse_blocks(raster,extra_raster)
% Collapse raster data (x, y, tr_idx) from different blocks.

for ir = 1:numel(raster)
   
    this_raster = raster(ir);
    
    rs_match=[];
    rs_match = find( [extra_raster.HP] == this_raster.HP ...
        & [extra_raster.LP] == this_raster.LP ...
        & [extra_raster.AMrate] == this_raster.AMrate ...
        & [extra_raster.dB] == this_raster.dB ...
        & [extra_raster.AMdepth] == this_raster.AMdepth ...
        & strcmp({extra_raster.behaving},this_raster.behaving) ...
        & strcmp([extra_raster.jitter],this_raster.jitter) ...
        );
    
    if numel(rs_match)>1, keyboard, end
    
    if isempty(rs_match)
        collapsed_raster(ir) = this_raster;
        continue
    end
    
    collapsed_raster(ir) = this_raster;
    % Add spike time info
    collapsed_raster(ir).x = [this_raster.x extra_raster(rs_match).x];
    collapsed_raster(ir).y = [this_raster.y extra_raster(rs_match).y + max(this_raster.y)*ones(size(extra_raster(rs_match).y)) ];
    % Add trial index info
    tr_idx = nan(50,2);
    tr_idx(1:length(this_raster.tr_idx),1) = this_raster.tr_idx;
    tr_idx(1:length(extra_raster(rs_match).tr_idx),2) = extra_raster(rs_match).tr_idx;
    collapsed_raster(ir).tr_idx = tr_idx(sum(~isnan(tr_idx),2)>0,:);
        
end

if numel(collapsed_raster) ~= numel(raster)
    keyboard
end

end
