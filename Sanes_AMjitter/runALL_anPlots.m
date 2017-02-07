function runALL_anPlots(subject,session,channel,clu)
% Plot all response measures for a designated unit.


% close all
% 
% METRICS = {'FR' 'FF' 'FFpd-avg' 'VS' 'RS' 'standardFR'};
% for im = 1:numel(METRICS)
%     ap_barplot_indJitter(subject,session,channel,clu,METRICS{im})
% end

METRICS = {'FF' 'VS' 'RS'};
for im = 1:numel(METRICS)
    ap_heatplot_indJitter(subject,session,channel,clu,METRICS{im})
end

% METRICS = {'FR'};
% for im = 1:numel(METRICS)
%     ap_neurometric_indDepth(subject,session,channel,clu,METRICS{im},0)
% end


end
