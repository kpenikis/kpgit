function [meanR_IR,meanR_pdc] = corr_comp_mph(Data_IR,Data_pdc)


for ii = 1:size(Data_pdc,1)
    
    % Get MPH to save to data table
    raster_x = []; raster_y = [];
    for it = 1:size(Data_pdc(ii,:).raster{:},1)
        raster_x = [raster_x find(Data_pdc(ii,:).raster{:}(it,:))];
        raster_y = [raster_y repmat(it,1,numel(find(Data_pdc(ii,:).raster{:}(it,:))))];
    end
    
    this_mph{ii,1} = hist(raster_x,linspace(0,1000/Data_pdc(ii,:).AMrate,52));
    
end

Data_pdc.MPH = this_mph;
clear this_mph

for ii = 1:size(Data_IR,1)
    
    % Get MPH to save to data table
    raster_x = []; raster_y = [];
    for it = 1:size(Data_IR(ii,:).raster{:},1)
        raster_x = [raster_x find(Data_IR(ii,:).raster{:}(it,:))];
        raster_y = [raster_y repmat(it,1,numel(find(Data_IR(ii,:).raster{:}(it,:))))];
    end
    
    this_mph{ii,1} = hist(raster_x,linspace(0,1000/Data_IR(ii,:).AMrate,52));
    
end

Data_IR.MPH = this_mph;


% Compute correlations across each pair of MPHs
[Rmat_IR,pmat_IR]  = corrcoef(vertcat(Data_IR.MPH{:})');
[Rmat_pdc,pmat_pdc] = corrcoef(vertcat(Data_pdc.MPH{:})');

% Find the average of these pairwise correlations
meanR_IR  = mean(Rmat_IR(Rmat_IR~=1));
meanR_pdc = mean(Rmat_pdc(Rmat_pdc~=1));



end






