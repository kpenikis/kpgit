
function stim = collapse_blocks(raster)
% Input only data to plot in single axes. (ie, already have separated
% different noise bands and SPLs).
% Collapse raster data from different blocks, keeping separate data based
% on jitter file ID and behavioral state.

drinking = strcmp({raster.behaving},'D');
behaving = strcmp({raster.behaving},'A');

[unq,~,iu] = unique([[raster.fileIDs]; drinking; behaving]','rows');
for iiu = 1:max(iu)
    
    rs = find([raster.fileIDs]==unq(iiu,1) & drinking==unq(iiu,2) & behaving==unq(iiu,3));
    
    x=[]; y=0; bk=[]; 
    for irs = rs
        x = [x raster(irs).x];
        y = [y raster(irs).y + max(y)*ones(size(raster(irs).y))];
        bk = [bk raster(irs).block];
    end
    y(1)=[];
    
    stim(iiu) = raster(rs(1));
    stim(iiu).block = bk;
    stim(iiu).x = x; 
    stim(iiu).y = y;
    
end
end
