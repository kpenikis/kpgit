function figname = create_figure(ij,var1,var2,dpcat,plottype,xlimits,ylimits,figsize,binsize)

global Figs 


%% ALL JITTERS TOGETHER

if ~isfield(Figs,dpcat)
    
    Figs.(dpcat) = dpcat;
    eval(sprintf('%s = figure;',Figs.(dpcat)))
    set(gcf,'Position',figsize,'Nextplot','add');
    
    switch plottype(1:4)
        case 'comp'
            title(sprintf('%s %s during standard period\n%s units',var1,plottype,dpcat))
            xlabel(sprintf('%s periodic',var1))
            ylabel(sprintf('%s jitter',var1))
            plot([-1 100],[-1 100],'k','LineWidth',0.5)
            set(gca,'xlim',xlimits,'ylim',ylimits)
        case 'regr'
            title(sprintf('%s as a function of %s\n%s units',var1,var2,dpcat))
            xlabel(sprintf('%s',var2))
            ylabel(sprintf('%s',var1))
            plot([-1 100],[0 0],'k','LineWidth',0.5)
            set(gca,'xlim',xlimits,'ylim',ylimits)
        case 'psth'
            xlim([1 250/binsize])
            title(sprintf('%s - during standard 4 Hz period\n%s units',plottype,dpcat))
            xlabel('time (ms) in standard pd')
            ylabel('Spikes/sec')
    end
end



%% EACH JITTER SEPARATELY

ijf = ij;
ijf(strfind(ijf,'-'))=[];

figname = [dpcat '_' ijf];

if ~isfield(Figs,figname)
    
    Figs.(figname) = figname;
    eval(sprintf('%s = figure;',Figs.(figname)))
    set(gcf,'Position',figsize,'Nextplot','add');
    
    switch plottype(1:4)
        case 'comp'
            title(sprintf('%s\n%s %s during standard period\n%s units',ij,var1,plottype,dpcat))
            xlabel(sprintf('%s periodic',var1))
            ylabel(sprintf('%s jitter',var1))
            plot([-1 100],[-1 100],'k','LineWidth',0.5)
            set(gca,'xlim',xlimits,'ylim',ylimits)
        case 'regr'
            title(sprintf('%s\n%s as a function of %s\n%s units',ij,var1,var2,dpcat))
            xlabel(sprintf('%s',var2))
            ylabel(sprintf('%s',var1))
            plot([-1 100],[0 0],'k','LineWidth',0.5)
            set(gca,'xlim',xlimits,'ylim',ylimits)
        case 'psth'
            xlim([1 250/binsize])
            title(sprintf('%s\n%s - during standard 4 Hz period\n%s units',ij,plottype,dpcat))
            xlabel('time (ms) in standard pd')
            ylabel('Spikes/sec')
            
    end
    
end


end



