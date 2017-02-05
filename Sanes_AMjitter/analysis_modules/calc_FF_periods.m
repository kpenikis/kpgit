
function [FFmean] = calc_FF_periods(stim,subject)

% Set stimulus file directory and load rate vectors
blocks = stim.block;
stimdir = fullfile('/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData',subject,sprintf('Block-%i_Stim',blocks(1)));
rV = load(fullfile(stimdir,stim(1).stimfn));

% Find minimum number of trials for each stimulus
ntc=nan(numel(stim),1);
for is = 1:numel(stim)
    ntc(is) = max(stim(is).y);
end
min_nt = min(ntc);
% min_nt = 5;

iterations = 200;  rng('shuffle');

% Set up empty vectors and get FF data
FF = nan( numel(stim), length(rV.buffer)-1, iterations );

for is = 1:numel(stim)
    
    data = stim(is);
    
    % Get vectors of rates for this stimulus
    rateVec = load(fullfile(stimdir,data.stimfn));
    rateVec = rateVec.buffer;
    tVec = round(data.AMonset + cumsum([0 0.75*(1000/rateVec(2)) 1000./rateVec(3:end)]));
    
for ii = 1:iterations
    
    kts = randperm(max(data.y));
    trs = kts(1:min_nt);
    sp_hist = nan(min_nt,length(tVec)-1);
    
    ir = 0;
    for it = trs
        ir=ir+1;
        
        sp = data.x(data.y==it);
        sp_hist(ir,:) = histcounts(sp,tVec);
        
    end
    
    % Calculate means and variances for time bins
    tr_var  = var(sp_hist,1);
    tr_mean = mean(sp_hist,1);
    
    % Calculate Fano Factor for this stimulus
    try
        tr_var(tr_mean==0) = nan; tr_mean(tr_mean==0) = nan;
        FF(is,:,ii) = tr_var(tr_mean~=0) ./ tr_mean(tr_mean~=0);
    catch
        keyboard
    end

end
end

FFmean = mean(FF,3);
FFstd  = std(FF,1,3)/sqrt(iterations);

end