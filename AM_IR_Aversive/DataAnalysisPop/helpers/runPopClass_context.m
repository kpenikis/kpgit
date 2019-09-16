% Called by PCl_MPcontext
% Goes through each MP comparison for each AM rate, for "theseUnits" 
% Run classifier on concatenated population data. 
% 

Pop_dPrimes = nan(4,5);
mintrs      = nan(1,5);
padTrs      = 300;

for irate = 1:5
    
    fprintf('  getting data for %i Hz...\n',AMrates(irate))
    
    PdcRaster_7_1 = nan(padTrs,0);
    PdcRaster_7_2 = nan(padTrs,0);
    PdcRaster_8_1 = nan(padTrs,0);
    PdcRaster_8_2 = nan(padTrs,0);
    
    IrrRaster_7_1 = nan(padTrs,0);
    IrrRaster_7_2 = nan(padTrs,0);
    IrrRaster_8_1 = nan(padTrs,0);
    IrrRaster_8_2 = nan(padTrs,0);
    
    
    % Step through Units and collect all spikes
    
    h = waitbar(0,sprintf('Concatenating raster data from all units for %i Hz MPs',AMrates(irate)));
    
    for iUn = theseUnits
        
        % Determine which Irr MPs were presented to this unit
        % now it's a field
        if iUn==0 || ~isfield(Data(iUn,irate).data,'IRid')
            fprintf('skipping iUn %i\n',iUn)
            continue
        end
        
        waitbar(iUn/numel(theseUnits))
        
        theseIrrMPs = [[Data(iUn,irate).data(:,1).IRid]; [Data(iUn,irate).data(:,1).iseq]]';
        
        for ic = 1:size(theseIrrMPs,1)
            
            % Skip if not enough trials (in either Pdc or Irr MP)
            if size(Data(iUn,irate).data(ic,1).raster,1)<trMin*2 || size(Data(iUn,irate).data(ic,2).raster,1)<trMin*2
                continue
            end
            
            % Skip if high d' and excluding best SU datapoints
            if exist('q','var') && Res(iUn,irate).L1o.dprime(ic,2)>q(irate)
                continue
            % Skip if LOW d' and keeping only best SU datapoints
            elseif exist('k','var') && Res(iUn,irate).L1o.dprime(ic,2)<k(irate)
                continue
            end
            
            
            % Concatenate Pdc raster data from this cell
            raster = nan(padTrs,ceil(1000/AMrates(irate)));
            raster(1:size(Data(iUn,irate).data(ic,1).raster,1),:) = Data(iUn,irate).data(ic,1).raster;
            
            eval(sprintf( 'PdcRaster_%i_%i = [PdcRaster_%i_%i raster];',theseIrrMPs(ic,1),theseIrrMPs(ic,2),theseIrrMPs(ic,1),theseIrrMPs(ic,2))) 
            
            
            % Concatenate Irr raster data from this cell
            raster = nan(padTrs,ceil(1000/AMrates(irate)));
            raster(1:size(Data(iUn,irate).data(ic,2).raster,1),:) = Data(iUn,irate).data(ic,2).raster;
            
            eval(sprintf( 'IrrRaster_%i_%i = [IrrRaster_%i_%i raster];',theseIrrMPs(ic,1),theseIrrMPs(ic,2),theseIrrMPs(ic,1),theseIrrMPs(ic,2))) 
            
        end %ic
        
    end %iUn
    
    close(h);
    
    % Trim NaNs
    PdcRaster_7_1 = PdcRaster_7_1(1:sum(~isnan(mean(PdcRaster_7_1,2))),:);
    PdcRaster_7_2 = PdcRaster_7_2(1:sum(~isnan(mean(PdcRaster_7_2,2))),:);
    PdcRaster_8_1 = PdcRaster_8_1(1:sum(~isnan(mean(PdcRaster_8_1,2))),:);
    PdcRaster_8_2 = PdcRaster_8_2(1:sum(~isnan(mean(PdcRaster_8_2,2))),:);
    IrrRaster_7_1 = IrrRaster_7_1(1:sum(~isnan(mean(IrrRaster_7_1,2))),:);
    IrrRaster_7_2 = IrrRaster_7_2(1:sum(~isnan(mean(IrrRaster_7_2,2))),:);
    IrrRaster_8_1 = IrrRaster_8_1(1:sum(~isnan(mean(IrrRaster_8_1,2))),:);
    IrrRaster_8_2 = IrrRaster_8_2(1:sum(~isnan(mean(IrrRaster_8_2,2))),:);
    
    % Get min n trials
    mintrs(irate) = min([size(PdcRaster_7_1,1) size(PdcRaster_8_1,1) size(IrrRaster_7_1,1) size(IrrRaster_8_1,1)]);
    
    % Put population data into struct
    clear PopMPH
    PopMPH = struct('Context',{'Pdc' 'IR'; 'Pdc' 'IR'; 'Pdc' 'IR'; 'Pdc' 'IR'},...
        'IRid',{7 7; 7 7; 8 8; 8 8},...
        'iseq',{1 1; 2 2; 1 1; 2 2},...
        'nTrs',{size(PdcRaster_7_1,1)  size(IrrRaster_7_1,1); ...
                size(PdcRaster_7_2,1)  size(IrrRaster_7_2,1); ...
                size(PdcRaster_8_1,1)  size(IrrRaster_8_1,1); ...
                size(PdcRaster_8_2,1)  size(IrrRaster_8_2,1) },...
        'raster',{PdcRaster_7_1 IrrRaster_7_1;...
                  PdcRaster_7_2 IrrRaster_7_2;...
                  PdcRaster_8_1 IrrRaster_8_1;...
                  PdcRaster_8_2 IrrRaster_8_2});
        %First draft, clip rasters to min N trials (later randomized which selected on each iteration)
    
    % Run classifier
    Pop_Res_L1o  = get_classifier_data( PopMPH, -1 );
    
    Pop_dPrimes(:,irate) = Pop_Res_L1o.dprime(:,2);
    
end %irate

save(fullfile(fn.processed,'MPcontext','Pop',['Pop_dPrimes_' append_str]),'Pop_dPrimes','-v7.3')


