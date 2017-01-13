function ap_plot_MPH(subject,session,channel,clu,raster)
% Plots a folded histogram of spiking during the AM portion of a stimulus.
% Note that periods are not all of equal length so n trials tapers off.

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

if nargin<5
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end
    
% Load new rate vector stimulus info
block = raster(1).block;
stimdir  = '/Users/kpenikis/Documents/SanesLab/Data/raw_data';
stfs = dir(fullfile(stimdir,subject,sprintf('Block-%i_Stim',block),'*.mat'));
stfns = {stfs.name};

% Get vectors of rates for stimuli of this block
for iif = 1:numel(stfs)
    rateVecs(iif) = load(fullfile(stfs(iif).folder,stfns{iif}));
end

for ks = 1:numel(raster)
    
    data = raster(ks);
    
    if block ~= data.block
        
        % Load new rate vector stimulus info
        block = data.block;
        stfs = dir(fullfile(stimdir,subject,sprintf('Block-%i_Stim',block),'*.mat'));
        stfns = {stfs.name};
        
        % Get vectors of rates for stimuli of this block
        for iif = 1:numel(stfs)
            rateVecs(iif) = load(fullfile(stfs(iif).folder,stfns{iif}));
        end
        
    end
    
    % Get rate vector for current stimulus
    rV = rateVecs(strcmp(stfns,data.stimfn)).buffer;
    
    % Define time vector for polar histo, boundaries between periods where
    % the rate changes
    periods = 1./rV(3:end-1)*1000;
    pd_range = [min(periods) max(periods)];
    t1 = 0.75*(1/rV(2)*1000) + data.AMonset; %ms
    pVec = round(cumsum([t1 periods]));
    
    % Set up for creating histo
    binsize = 10;
    tVec = binsize/2 : binsize : ceil(pd_range(2)/binsize)*binsize;
    nt = max(data.y);
    sp_histo = nan( numel(pVec)-1, ceil(pd_range(2)/binsize) );
    sp_histo_avg = nan( nt, ceil(pd_range(2)/binsize) );
    
    % Begin figure
    hF(ks) = figure;
    fill( [pd_range(1) pd_range(1) pd_range(2) pd_range(2)], [0 300 300 0],...
        [0.95 0.95 0.95], 'EdgeColor','none' )
    hold on
    
    % Cycle through trials and get spiketimes
    for it = 1:nt
        sp=[]; sp=data.x(data.y==it);
        
        % Cycle through each period of AM
        for ip = 1:numel(pVec)-1
            
            % Define edges for histogram (according to nbins set above)
            tp = pVec(ip):binsize:pVec(ip+1);
            hc=[]; hc = histcounts(sp(sp>pVec(ip)&sp<(pVec(ip+1))),tp)/(binsize/1000);
            sp_histo(ip,1:length(hc)) = hc;
            
        end
                
        % Average "MPH" for this trial
        sp_histo_avg(it,:) = sum(sp_histo,1,'omitnan')/ip;
        plot(tVec, sp_histo_avg(it,:))
        
    end
    
    % Plot average for the stimulus
    plot(tVec,sum(sp_histo_avg,1,'omitnan')/nt,'k','LineWidth',2)
    ylim([0 max(max(sp_histo_avg))])
    xlim([0 pd_range(2)])
    ylabel('Firing rate (sp/s)')
    xlabel('Time within wrapped AM periods (ms)')
    title(raster(ks).stim_str)
    hold off
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    savename = sprintf('%s_%s_MPH_ch%i_clu%i_AM4Hz_jitter%s_blk%i',subject,session,channel,clu,data.jitter,block);
    print(hF(ks),'-depsc',fullfile(an_dir,savename))
    
end

end

