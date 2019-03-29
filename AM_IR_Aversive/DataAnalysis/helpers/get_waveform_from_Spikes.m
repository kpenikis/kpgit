
% Create intepolated version of waveform
spikes = Spikes.sorted(channel);
waveform = median(spikes.waveforms(spikes.assigns==clu,:),1);
waveform = waveform./abs(min(waveform));
m = 40;
x = (1:length(waveform)) /Info.fs*1000;
q = linspace(min(x),max(x),length(x)*m);
waveform_interp = interp1(x,waveform,q,'spline');


% Now, if this was a merged unit, load the other waveform and compare
% to the first
if numel(sess_tmp)>1
    
    session2 = sess_tmp{2};
    clu2     = UnitData(iUn).Clu(2);
    filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session2); load(fullfile(fn.processed,subject,filename));
    
    % Create intepolated version of waveform
    spikes = Spikes.sorted(channel);
    waveform2 = median(spikes.waveforms(spikes.assigns==clu2,:),1);
    waveform2 = waveform2./abs(min(waveform2));
    m = 40;
    x = (1:length(waveform2)) /Info.fs*1000;
    q = linspace(min(x),max(x),length(x)*m);
    waveform_interp2 = interp1(x,waveform2,q,'spline');
    
    %         hf_tmp2 = figure; hold on
    %         plot(q,waveform_interp,'r','LineWidth',2)
    %         plot(q,waveform_interp2,'k','LineWidth',2)
    %
    %         keyboard
    %         close(hf_tmp2);
    
    waveform        = mean([waveform; waveform2],1);
    waveform_interp = mean([waveform_interp; waveform_interp2],1);
end


