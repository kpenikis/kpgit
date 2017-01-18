
function [FFmean,FFpctl] = calc_FF(stim,binsize)

% Find minimum number of trials for each stimulus
ntc=nan(numel(stim),1);
for is = 1:numel(stim)
    ntc(is) = max(stim(is).y);
end
min_nt = min(ntc);
% min_nt = 8;

iterations = 200;  rng('shuffle');

% Set up empty vectors and get FF data
FF     = nan(numel(stim),iterations);

for is = 1:numel(stim)
    
    data = stim(is);
    
    if binsize==0
        tVec = [data.AMonset data.stimDur];
    else
        tVec = data.AMonset : binsize : data.stimDur;
    end
    
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
        FF(is,ii) = mean(tr_var(tr_mean~=0) ./ tr_mean(tr_mean~=0));

    end
end


FFmean = mean(FF,2);
FFster = [std(FF,1,2)/sqrt(iterations) std(FF,1,2)/sqrt(iterations)] ./2;
FFpctl = [FFmean-prctile(FF,25,2) prctile(FF,75,2)-FFmean];

end

