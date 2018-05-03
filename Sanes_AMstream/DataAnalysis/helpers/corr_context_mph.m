function R_context = corr_context_mph(Data_IR,Data_pdc,thisAMrate)


[~,row_idx] = sort(Data_IR.Starttime);
R_context   = nan(numel(row_idx),2);

for irow = 1:numel(row_idx)
    
    %~~~~~~~~~~
    % IR first
    
    % Get raster data
    raster_x = [];
    for it = 1:size(Data_IR(row_idx(irow),:).raster{:},1)
        raster_x = [raster_x find(Data_IR(row_idx(irow),:).raster{:}(it,:))];
    end
    % Get MPH
    MPH_IR = hist(raster_x,linspace(0,1000/thisAMrate,52));
    
    
    %~~~~~~~~~~~~~~
    % Now Periodic
    
    [~,this_pdc] = min(abs(Data_pdc.Starttime - Data_IR(row_idx(irow),:).Starttime));
    
    % Get raster data
    raster_x = [];
    for it = 1:size(Data_pdc(this_pdc,:).raster{:},1)
        raster_x = [raster_x find(Data_pdc(this_pdc,:).raster{:}(it,:))];
    end
    % Get MPH
    MPH_pdc = hist(raster_x,linspace(0,1000/thisAMrate,52));
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Calculate correlation coefficient
    [r,p] = corrcoef(MPH_IR,MPH_pdc);
    R_context(irow,1) = r(2,1);
    R_context(irow,2) = p(2,1);
    
    
end




end


