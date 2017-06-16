
function [FRvec,FRstd,FRtr] = calc_standardFR(stim,subject)
% now vertical output vector

global fn

% Set stimulus file directory and load rate vectors
blocks = stim.block;
stimdir = fullfile(fn.raw,subject,sprintf('Block-%i_Stim',blocks(1)));

% Get FR and std of FR for error bars
FRvec = nan(numel(stim),1); FRstd = nan(numel(stim),2); FRtr=struct();
for is = 1:numel(stim)
    
    data = stim(is);
    
    % if FR is too low, set data output to nans
    if isempty(data.x) || (numel(data.x)/max(data.y) / (data.stimDur/1000)) < 5
        disp('skipping datapoint with too few spikes')
        continue  
    end
    
    % Get vectors of rates for this stimulus
    rateVec = load(fullfile(stimdir,data.stimfn));
    rateVec = rateVec.buffer(2:end);
    pd_standard = 4;%ceil((length(rateVec))/2);
    
    t0 = ceil(data.AMonset + 0.75*1000/rateVec(1) + sum(1000./rateVec(2:pd_standard-1)));
    tVec = [t0 t0+1000/rateVec(pd_standard)];
    
    if diff(tVec)~=250
        keyboard
    end
    
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
