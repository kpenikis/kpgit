function ap_plot_polarHisto(subject,session,channel,clu,raster)
% Plots a polar histogram of spiking during the AM portion of a stimulus. 

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
    hF(ks) = figure;
    
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
    t1 = 0.75*(1/rV(2)*1000) + data.AMonset; %ms
    tVec = round(cumsum([t1 periods]));
    
    % Set up for creating histo
    nbins = 360/2;
    theta = [2*pi/nbins : 2*pi/nbins : 2*pi] - (2*pi/nbins/2);
    nt = max(data.y);
    sp_phase_counts = nan(numel(tVec)-1,nbins);
    sp_phase_avg = nan(nt,nbins);
    
    % Cycle through trials and get spiketimes
    for it = 1:nt
        sp=[]; sp=data.x(data.y==it);
        
        % Cycle through each period of AM
        for ip = 1:numel(tVec)-1
            
            % Define edges for histogram (according to nbins set above)
            tp = round(linspace(tVec(ip),tVec(ip+1),nbins+1));
            sp_phase_counts(ip,:) = histcounts(sp(sp>=tVec(ip)&sp<(tVec(ip+1))),tp);
            
        end
        
        % Average "MPH" for this trial
        sp_phase_avg(it,:) = sum(sp_phase_counts,1);
        
        polarplot(theta,sum(sp_phase_counts,1))
        hold on
        
    end
    
    % Get VS for this stim
    sp_avg = sum(sp_phase_avg,1)/nt;
    keyboard
    VS = theta .* sp_avg;
    meanx = atan2(mean(sin(VS)),mean(cos(VS)))
    VS = mean(VS)
    
    % Plot average for the stimulus
    polarplot(theta,sum(sp_phase_avg,1)/nt,'k','LineWidth',2)
    title(raster(ks).stim_str)
    polarplot([meanx meanx],[0 4*max(sp_avg)],'r-','LineWidth',3)
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    savename = sprintf('%s_%s_polarHisto_deg_ch%i_clu%i_AM4Hz_jitter%s_blk%i',subject,session,channel,clu,data.jitter,block);
    print(hF(ks),'-depsc',fullfile(an_dir,savename))
    
end




end

