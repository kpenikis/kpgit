
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/MPHclassifier/ClassData.mat')
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/Units.mat')

for irate = 2
    
    for iUn = 1:size(Data,1)
        
        idx = Data(iUn,irate).Res_L1o.dprime(Data(iUn,irate).Res_L1o.dprime(:,2) > 2, 1);
        
        if isempty(idx) || mean(mean(Data(iUn,irate).data(1).raster)) > 1 
            continue
        end
                
        for ii = idx
            
            raster_Irr = Data(iUn,irate).data(ii).raster;
            raster_Pdc = Data(iUn,irate).data(1).raster;
            
            hf=figure;
            plot(mean(raster_Pdc,1),'k')
            hold on
            plot(mean(raster_Irr,1),'b')
            
            keyboard
            
        end
        
    end
    
end


