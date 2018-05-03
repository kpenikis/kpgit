function UnitPlotMetrics(subject,session,channel,Clus)
%%
% Load data files
fn = set_paths_directories(subject,session);
fprintf('Session %s:  loading data...\n',session)
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));



% Find clus to plot
spikes = Spikes.sorted(channel);

if nargin<4
    Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
end

for clu = Clus'
    
    figure; 
    subplot(2,2,1)
    plot([0 38], [0 0], 'k','LineWidth',0.25)
    hold on
    fill([1:size(spikes.waveforms,2) fliplr(1:size(spikes.waveforms,2))],[mean(spikes.waveforms(spikes.assigns==clu',:),1)+std(spikes.waveforms(spikes.assigns==clu',:),1) fliplr(mean(spikes.waveforms(spikes.assigns==clu',:),1)-std(spikes.waveforms(spikes.assigns==clu',:),1))],...
        'k','FaceAlpha',0.3,'EdgeColor','none')
    plot(mean(spikes.waveforms(spikes.assigns==clu',:),1),'k','LineWidth',1)
    ylim([-6 2].*10^-4)
    xlim([1 37])
    
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    timebins = 0:51;
    subplot(2,2,2)
    halfhist = hist(diff(spiketimes),timebins);
    bar([-fliplr(timebins(2:end-1)) timebins(1:end-1)],[fliplr(halfhist(2:end-1)) halfhist(1:end-1)],'k')
    xlim([-timebins(end-1) timebins(end-1)])
    
    keyboard
    print_eps_kp(gcf,fullfile(fn.processed,'UnitPlots',sprintf('WaveformAutocorr_%s_ch%i',session,channel,clu)),1)
    
end


% For PC plot, launch UMS gui and run the next 2 lines to save
pp_launch_manual_sort(subject,session,channel);

%%

axis square
print_eps_kp(gcf,fullfile(fn.processed, 'UnitPlots', sprintf('PCplot_%s_ch%i',session,channel)),1)


%%



end %function
