function FanoFs = FanoFactorDistributions(plot_flag)
% 
% Add: var (mean) n spikes, by stimulus (log space) 
% First: get distribution of RTs 
%

global AMrates 

% alfaVS    = 0.001;
% trMax     = 20;

Stimuli   = 1:8;
Duration  = 1000;
skipOnset = 0;
getMPH    = 0;
        
% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
twothirds   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)];
widescreen  = [1 scrsz(4) scrsz(3) scrsz(4)/3*2];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];



%%

if plot_flag
    hfsp=figure;
    hold on
end

FanoFs  = zeros(numel(UnitData),numel(Stimuli));

for iUn = 1:numel(UnitData)
    
    get_trial_data_posthoc
    
    FanoFs(iUn,:) = [StimSpikeData(:,2)./StimSpikeData(:,1)]';
    
    if plot_flag
        plot(StimSpikeData(1,1),StimSpikeData(1,2),'ok')
    end
    
end %iUn

% Save FF data
save(fullfile(fn.figs,'FF',sprintf('FF_Stim%i-%i_%i-%ims',min(Stimuli),max(Stimuli),skipOnset,skipOnset+Duration)),'FanoFs','-v7.3')


if plot_flag
    % Finish n spikes plot
    axis square
    plot(10.^[-2 3],10.^[-2 3],'--k')
    ylim(10.^[-2 3])
    xlim(10.^[-2 3])
    title(sprintf('Stim %i-%i: %i-%i ms',min(Stimuli),max(Stimuli),skipOnset,skipOnset+Duration))
    xlabel('Mean n spikes')
    ylabel('Var n spikes')
    set(gca,'xscale','log','yscale','log',...
        'xtick',[0.01 0.1 1 10 100 1000],'xticklabel',[0.01 0.1 1 10 100 1000],...
        'ytick',[0.01 0.1 1 10 100 1000],'yticklabel',[0.01 0.1 1 10 100 1000])
    set(gca,'tickdir','out','ticklength',[0.04 0.04],'Color','none')
    
    print_eps_kp(hfsp,fullfile(fn.figs,'FF',sprintf('Nspk_Stim%i-%i_%i-%ims',min(Stimuli),max(Stimuli),skipOnset,skipOnset+Duration)))
end


%% Plot results

if plot_flag
    
    if size(FanoFs,2)==9
        FanoFs        = [ones(size(FanoFs,1),1) FanoFs];
        FanoFs(:,1)   = FanoFs(:,end);
        FanoFs(:,end) = [];
        xaxislabels = ['Spont.' Info.stim_ID_key' ];
    else 
        xaxislabels = Info.stim_ID_key(Stimuli)';
    end
    
    hf1 = figure;
    set(gcf,'Position',fullscreen)
    hold on
    
    FFmeds = nan(1,5);
    for ist = 1:size(FanoFs,2)
        
        FFmeds(ist) = median(FanoFs(~isnan(FanoFs(:,ist)),ist));
        
        % Manually make boxplot
        q5  = quantile(FanoFs(:,ist),0.05);
        q25 = quantile(FanoFs(:,ist),0.25);
        q75 = quantile(FanoFs(:,ist),0.75);
        q95 = quantile(FanoFs(:,ist),0.95);
        
        plot(ist*[1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
        fill(ist+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
        
    end
    
    plotSpread(FanoFs,'showMM',0,'distributionColors','k')
    plot(1:size(FanoFs,2),FFmeds,'.-','MarkerSize',40,'LineWidth',2,'Color',[1 0.8 0])
    
    set(gca,'xtick',1:numel(xaxislabels),'xticklabel',xaxislabels)
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    xlim([0 10])
    ylim([0 10])
    xlabel('AM rate (Hz)')
    ylabel('FF')
    title(sprintf('FF: %i-%i ms',skipOnset,skipOnset+Duration))
    
    axis fill
    box on
    
    print_eps_kp(hf1,fullfile(fn.figs,'FF',sprintf('FF_Stim%i-%i_%i-%ims',min(Stimuli),max(Stimuli),skipOnset,skipOnset+Duration)))
    
end


end





