function VSdata_pdc = get_VS_periodic(BinarySpikeData,AMrates)
% VSdata_pdc = get_VS_periodic(BinarySpikeData,AMrates)
%  Called by ContextMTFs. Only for periodic stimuli.
%  Output is a structure. 
%  Input is an 8x8 cell with rasters in sparse matrix format.
% 
% KP, 2018-04
%

VSdata_pdc = struct;

for ist = 2:6
    
    % Get minimum number of trials
    Ntrs = cell2mat( cellfun(@(x) size(x,1) , BinarySpikeData(ist,:),'UniformOutput',false));
    minNtrs = min(Ntrs(Ntrs~=0));
    
    period = 1000/AMrates(ist-1);
    
    for ipst = [3 6]
        
        sp01 = BinarySpikeData{ist,ipst};
        if isempty(sp01)
            VSdata_pdc(ist-1,ipst/3).VS = nan;
            VSdata_pdc(ist-1,ipst/3).RS = nan;
            VSdata_pdc(ist-1,ipst/3).RP = nan;
            continue
        end
        
        VS = nan(1,size(sp01,1));
        RS = nan(1,size(sp01,1));
        RP = nan(1,size(sp01,1));
        
%         if size(sp01,1)>minNtrs
                % bootstrap trials and average the results, to stabilize
%         else
            for it = 1:size(sp01,1)
                spktimes = find(sp01(it,:));
                [VS(it),RS(it),RP(it)] = vectorstrength(spktimes,period);
            end
%             if any(isnan(VS))
%                 keyboard
%             end
            VSdata_pdc(ist-1,ipst/3).VS = mean(VS);
            VSdata_pdc(ist-1,ipst/3).RS = mean(RS);
            VSdata_pdc(ist-1,ipst/3).RP = mean(RP);
%         end
        
        clear VS; clear RS; clear RP
        
    end %ipst
end %ist


% Check if unit is significantly synchronized, according to criterion from
% Malone 2007.
if ~any([VSdata_pdc(1:5,:).RP]<0.05)
%     keyboard
end


end

