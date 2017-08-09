function ProcessTuningSession(subject,sessionID,blocks)

pp_prepare_format( blocks, subject, sessionID )

pp_sort_tuning_session( subject, sessionID )

ap_plot_TC_AMrate( subject, sessionID, channel )
%change name

end
