
function [FRvec,FRstd,FRtr] = calc_FR(stim)
% [FRvec,FRstd,FRtr] = calc_FR(stim)
%  Get FR and std of FR, formatted for error bars.
%  data output is now vertical vector
%  KP 2017
%

FRvec = nan(numel(stim),1); FRstd = nan(numel(stim),2); FRtr=struct();

for ir = 1:numel(stim)
    
    % if FR is too low, set data output to nans
    if isempty(stim(ir).x) || (numel(stim(ir).x)/max(stim(ir).y) / (stim(ir).stimDur/1000)) < 5
        disp('skipping datapoint with too few spikes')
        continue 
    end
    
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