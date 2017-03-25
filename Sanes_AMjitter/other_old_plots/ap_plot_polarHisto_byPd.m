function ap_plot_polarHisto_byPd(subject,session,channel,clu,raster)
% Plots a polar histogram of spiking during the AM portion of a stimulus. 

% IF SAVING PDF FILES
figFontSize      = 16;
rasterMarkerSize = 10;
rasterLineWidth  = 0.5;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;


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
    
    if numel(raster(ks).tr_idx)<6
        continue
    end
    
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
    periods = 1./rV(3:end)*1000;
    t1 = 0.75*(1/rV(2)*1000) + data.AMonset; %ms
    tVec = round(cumsum([t1 periods]));
    
    % Set up for creating histo
    nbins = 25;
    theta = [2*pi/nbins : 2*pi/nbins : 2*pi] - (2*pi/nbins/2);
    nt = max(data.y);
    sp_phase_avg = zeros(numel(tVec)-1,nbins);
    
    pdcolors = jet(numel(periods));
    
    % Cycle through each period of AM 
    for ip = 1:numel(tVec)-1
        
        % Define edges for histogram (according to nbins set above)
        tp = round(linspace(tVec(ip),tVec(ip+1),nbins+1));
                
        % Cycle through trials and get spiketimes
        sp_trial_counts = zeros(nt,nbins);
        for it = 1:nt
            
            sp=[]; sp=data.x(data.y==it);
            
            sp_trial_counts(it,:) = histcounts(sp(sp>tVec(ip)&sp<(tVec(ip+1))),tp);
            
        end
        
        % Average "MPH" for this trial
        sp_phase_avg(ip,:) = mean(sp_trial_counts,1);
        
        polarplot(theta,sp_phase_avg(ip,:) ./ ((tVec(ip+1)-tVec(ip))/1000/nbins),'Color',pdcolors(ip,:),'LineWidth',1.5)
        hold on
        
        legstr{ip} = ['pd ' num2str(ip)];
        
    end
    
    % Plot average for the stimulus
    try
    polarplot(theta,mean(sp_phase_avg./((tVec(ip+1)-tVec(ip))/1000/nbins),1),'k','LineWidth',3)
    title(data.stim_str)
    rlim([0 175])
    ax=gca;
    ax.ThetaZeroLocation = 'bottom';
    ax.ThetaDir = 'clockwise';
    legend(legstr,'mean')
    catch
        keyboard
    end
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    savename = sprintf('%s_%s_polarHisto-Pd_ch%i_clu%i_AM%iHz_jitter%s_%idpth_%idB_%i-%i_%s_blk%i',subject,session,channel,clu,...
        data.AMrate, data.jitter, round(data.AMdepth*100),...
        data.dB, data.HP, data.LP, data.behaving, data.block);
%     print(hF(ks),'-depsc',fullfile(an_dir,savename))
    set(gcf,'PaperOrientation','landscape');
    print(hF(ks),'-dpdf',fullfile(savedir,subject,'^an_plots',session,savename),'-bestfit')
    
end




end

