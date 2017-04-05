
function [FRvec,FRstd,FRtr] = calc_FR(stim)

% Get FR and std of FR for error bars
FRvec = nan(1,numel(stim)); FRstd = nan(numel(stim),2); FRtr=struct();
for ir = 1:numel(stim)
    x = stim(ir).x;
    y = stim(ir).y;
    FR_t = nan(1,max(y));
    for it = 1:max(y)
        FR_t(it) = numel(x(x(y==it) > stim(ir).AMonset & x(y==it) < stim(ir).stimDur)) / ((stim(ir).stimDur-stim(ir).AMonset)/1000);
    end
    FRtr(ir).tr = FR_t;
    FRvec(ir) = mean(FR_t,'omitnan');
    FRstd(ir,:) = [std(FR_t,'omitnan') std(FR_t,'omitnan')]./2;
end

end