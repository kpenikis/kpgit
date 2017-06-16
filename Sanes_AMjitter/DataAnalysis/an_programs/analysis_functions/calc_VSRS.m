function [VSmean,VS,RSmean,RS] = calc_VSRS(raster,subject,trials)

global fn

% Set stimulus file directory
blocks = raster.block;
stimdir = fullfile(fn.raw,subject,sprintf('Block-%i_Stim',blocks(1)));
rV = load(fullfile(stimdir,raster(1).stimfn));

% Set up empty vectors and get VS data
VS = nan(numel(raster),length(rV.buffer)-2);
RS = nan(numel(raster),length(rV.buffer)-2);
for is = 1:numel(raster)
    
    data = raster(is);
    
    % if FR is too low, set data output to nans
    if isempty(data.x) || (numel(data.x)/max(data.y) / (data.stimDur/1000)) < 5
        disp('skipping datapoint with too few spikes')
        continue  %VS already full of nans
    end
    
    % Get vectors of rates for this stimulus
    rateVec = load(fullfile(stimdir,data.stimfn));
    rateVec = rateVec.buffer;
    tVec = round(data.AMonset + cumsum([0.75*(1000/rateVec(2)) 1000./rateVec(3:end)]));
    
    % collect spikes across trials for each period
    if nargin<3
        trs = unique([data.y]);
    else
        trs = trials;
    end
    for ipd = 1:numel(tVec)-1
        AMpdms = (tVec(ipd+1)-tVec(ipd));
        
        Spks=[];
        for it = trs
            sp=[];
            sp = data.x(data.y==it);
            Spks = [Spks ( sp( sp>=tVec(ipd) & sp<=(tVec(ipd+1)) ) -tVec(ipd)) + ((it-1)*AMpdms)  ];
        end
        
        % Check for extremely low FR
        if numel(Spks)<numel(trs)
            continue
        end
        
        % for this period, call -vectorstrength-
        [VS(is,ipd),RS(is,ipd)] = vectorstrength(Spks,AMpdms,max(trs));
        
    end    
    
end

VSmean = mean(VS,2);
RSmean = mean(RS,2);



end