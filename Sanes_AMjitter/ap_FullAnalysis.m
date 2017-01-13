
function ap_FullAnalysis(subject,session,channel,clu)

% pp_plot_rasters_wStim( subject, session, channel, clu)
% pause(1)
% close all

ap_plot_MPH( subject, session, channel, clu)
pause(1)
close all

ap_plot_polarHisto_v1( subject, session, channel,clu)
pause(1)
close all

ap_plot_FRresp( subject, session, channel, clu)
ap_plot_FF( subject, session, channel, clu)


end

