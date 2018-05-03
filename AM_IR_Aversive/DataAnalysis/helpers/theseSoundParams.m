function [dBSPL,LP] = theseSoundParams(TrialData)
% [dBSPL,LP] = theseSoundParams(TrialData)
% Called by ap_plot_rasters, etc
% KP, 2018-03
%

dBSPL = unique(TrialData.SPL);
rm_i=[];
for ii = 1:numel(dBSPL)
    if sum(TrialData.SPL==dBSPL(ii))  < 300
        rm_i = [rm_i ii];
    end
end
dBSPL(rm_i) = [];

% Get unique noisebands (based on LP)
LP = unique(TrialData.LP);
rm_i=[];
for ii = 1:numel(LP)
    if sum(TrialData.LP==LP(ii))  < 300
        rm_i = [rm_i ii];
    end
end
LP(rm_i) = [];

end