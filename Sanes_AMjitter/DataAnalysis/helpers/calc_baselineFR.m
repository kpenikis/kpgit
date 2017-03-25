function baseFR = calc_baselineFR(raster)

if isempty(raster), baseFR=nan; return, end

nsp    =   sum([raster.x]<0);
t_s    =   (raster(1).window_ms(1)/-1000);
ntr    =   sum(cellfun(@(x) ( numel(x) ), {raster.tr_idx}));

baseFR = nsp / ntr / t_s;

end

    
    