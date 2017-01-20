
function [FRvec,FRstd,FRtr] = calc_standardFR(stim,subject)

% Set stimulus file directory and load rate vectors
blocks = stim.block;
stimdir = fullfile('/Users/kpenikis/Documents/SanesLab/Data/raw_data',subject,sprintf('Block-%i_Stim',blocks(1)));

% Get FR and std of FR for error bars
FRvec = nan(1,numel(stim)); FRstd = nan(numel(stim),2); FRtr=struct();
for is = 1:numel(stim)
    
    data = stim(is);
    
    % Get vectors of rates for this stimulus
    rateVec = load(fullfile(stimdir,data.stimfn));
    rateVec = rateVec.buffer;
    pd_standard = ceil((length(rateVec)-1)/2);
    t0 = ceil(data.AMonset + 0.75*1000/rateVec(3) + sum(1000./rateVec(3:pd_standard-1)));
    tVec = [t0 t0+1000/rateVec(pd_standard)];
    
    FR_t = [];
    trs = 1:max(data.y);
    for it = trs
        sp = data.x(data.y==it);
        sp_hist(it) = histcounts(sp,tVec);
    end
    
    FR_t = sp_hist*rateVec(pd_standard);

    FRtr(is).tr = FR_t;
    FRvec(is) = mean(FR_t,'omitnan');
    FRstd(is,:) = [std(FR_t,'omitnan') std(FR_t,'omitnan')]./2;
end

end
