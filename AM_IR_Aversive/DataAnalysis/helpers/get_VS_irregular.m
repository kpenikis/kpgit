function VSdata_OUT = get_VS_irregular(BinarySpikeData)
% get_VS_irregular = get_VS_periodic(BinarySpikeData,AMrates)
%  Called by ContextMTFs. For irregular stimuli.
%  Output is a structure. 
%  Input is an 8x8 cell with rasters in sparse matrix format.
% 
% KP, 2018-04
%

global fn

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

VSdata_OUT=struct;

for ist = 7:8
    
    switch ist
        case 7
            allPds = rateVec_AC;
        case 8
            allPds = rateVec_DB;
    end
    newPds = [1 cumsum(1000./allPds)];
    
    % Get minimum number of trials
    Ntrs = cell2mat( cellfun(@(x) size(x,1) , BinarySpikeData(ist,:),'UniformOutput',false));
    minNtrs = min(Ntrs(Ntrs~=0));
    
    for ipst = 2:6
        
        sp01 = BinarySpikeData{ist,ipst};
        if isempty(sp01)
            VSdata_OUT(ist-6,ipst-1).VS = nan;
            VSdata_OUT(ist-6,ipst-1).RS = nan;
            VSdata_OUT(ist-6,ipst-1).RP = nan;
            continue
        end
        
        VS = nan(1,numel(allPds));
        RS = nan(1,numel(allPds));
        RP = nan(1,numel(allPds));
        
        for ipd = 1:numel(allPds)
            
            period = 1000/allPds(ipd);
            
            spktimes=[];
            for it = 1:size(sp01,1)
                spktimes = [spktimes period*(it-1)+find(sp01(it,round(newPds(ipd):(newPds(ipd+1)-1))))];
            end
            
            [VS(ipd),RS(ipd),RP(ipd)] = vectorstrength(spktimes,period);
            
        end %ipd
        
        VSdata_OUT(ist-6,ipst-1).VS = VS;
        VSdata_OUT(ist-6,ipst-1).RS = RS;
        VSdata_OUT(ist-6,ipst-1).RP = RP;
        
    end %ipst
    
end %ist

end




