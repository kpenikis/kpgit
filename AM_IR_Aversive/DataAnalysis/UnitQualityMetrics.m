function UnitQualityMetrics(SUBJECT,SESSION)
%
%  UnitQualityMetrics(SUBJECT, SESSION )
%    Plots a raster and psth for each stimulus, for all SU from the session.
%    Excludes datapoints based on: min Ntrials, min FR.
%
%  KP, 2019-01
%

close all

%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];


%%

nfh   = 0;
nrws  = 5;
npans = nrws*2;

% for iUn = 1:numel(Clusters)
%     if mod(iUn,npans)==1
%         ip  = 0;
%     end
%     
%     mod(iUn,npans)+1
%     
% end
for iUn = 1:numel(Clusters)
    
    % Get spiketimes
    thisClu = Clusters(iUn); 
    maxChan = thisClu.maxChannel;
    
    %originally: seconds, not rounded
    spiketimes = round(thisClu.spikeTimes*1000)';
%     meanFR     = 1000* sum(Stream_Spikes(TrialData.onset(1):end)) / length(Stream_Spikes(TrialData.onset(1):end));
    baseFR = sum(spiketimes>TrialData.onset(1) & spiketimes<TrialData.offset(1)) / ((TrialData.offset(1)-TrialData.onset(1))/1000);
    
    % Plot unit summary
    
    if mod(iUn,npans)==1
        nfh = nfh+1;
        ip  = 0;
        hfClu(nfh) = figure;
        set(gcf,'Position',fullscreen)
        hold on
    end
    ip = ip+1;
    
    subplot(nrws,10, ip*5-4 + 0 );
    plot(thisClu.maxChTemp,'k','LineWidth',3)
    set(gca,'xlim',[1 length(thisClu.maxChTemp)],'xtick',[],'Color','none')
    ylim([-0.03 0.01])
    title(sprintf('spont %0.1f Hz',baseFR ))
    
    subplot(nrws,10, ip*5-4 + 1 );
    histogram([-diff(spiketimes) diff(spiketimes)],[-199.5:199.5],'EdgeColor','none','FaceColor','k','FaceAlpha',1);
    set(gca,'xlim',[-200 200],'Color','none')
    title(sprintf('%0.1f%% violations', 100*sum(diff(spiketimes)<3)/length(spiketimes) ))
    
    subplot(nrws,10, ip*5-4 + 2 );
    try
        histogram([thisClu.Amplitudes],50,'FaceColor','k','EdgeColor','none')
        xlim([0 50])
    end
    title(sprintf('%i:  ch %i  %i',thisClu.shank,maxChan,thisClu.clusterID ))
    
    try
        subplot(nrws,10, ip*5-4 + [3 4] );
        plot(spiketimes,[thisClu.Amplitudes],'.k')
        xlim([0 TrialData.offset(end)])
        title(sprintf('%i events',length(spiketimes)))
        ylim([0 50])
    catch
        keyboard
    end
        
    
end % iUn

%% Save figures

savedir = fullfile(fn.sessdata,'Rasters');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

for ihf = 1:nfh
    savename = sprintf('%s_%s_%i',SUBJECT,SESSION,ihf);
    print_eps_kp(hfClu(ihf),fullfile(savedir,savename))
%     set(hfClu(ihf),'PaperOrientation','landscape')
%     print(hfClu(ihf),fullfile(savedir,savename),'-dpdf','-bestfit')
end



end