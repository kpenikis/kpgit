function plotTransitionPSTH( IR_psth_wins, psth, stim, ContextString, TitleString, Durations, ymaxval, raster_x_1, raster_x_2, raster_y_1, raster_y_2 )

global t_win wincolors pstcolors figsize_psth 


% Make figure
hfr = figure;
set(hfr,'Position',figsize_psth,'NextPlot','add')
hold on
% Create subplots
clear hs ip
hs(1)=subplot(5,1,1);   box off; set(gca,'Color','none')
hs(2)=subplot(5,1,2:3); box off; set(gca,'Color','none')
hs(3)=subplot(5,1,4:5); box off; set(gca,'Color','none')
plot([0 0],[0 ymaxval],'k--')
hold on

% Plot the space between PSTHs
for iw = 2:numel(t_win)
    patch( -1+[t_win(iw-1):t_win(iw)-1 t_win(iw)-1:-1:t_win(iw-1)] ,...
        [IR_psth_wins(end,:,iw-1) fliplr(IR_psth_wins(1,:,iw-1))] , wincolors(iw-1,:),...
        'EdgeColor','none','FaceAlpha',0.65);
end

for iContext = 1:2
    
    subplot(hs(3)); hold on
    ip(iContext) = plot( -Durations(1):Durations(2), ...
        mean(psth(iContext,:,:),3,'omitnan') ,...
        'Color',pstcolors(iContext,:),'LineWidth',4);
    hold off
    legstr{iContext} = [ContextString{iContext} ', n=' num2str(max((eval(['raster_y_' num2str(iContext)]) )))];
    
    subplot(hs(1)); hold on
    plot(-Durations(1):Durations(2), mean(stim(iContext,:,:),3,'omitnan'),...
        'Color',pstcolors(iContext,:),'LineWidth',4)
    hold off
    
    add_y=0;
    if numel(eval(['raster_y_' num2str(iContext)]) )>0
        subplot(hs(2)); hold on
        plot(eval(['raster_x_' num2str(iContext)]),eval(['raster_y_' num2str(iContext)])+add_y,...
            '.','MarkerSize',15,'Color',pstcolors(iContext,:))
        add_y = add_y + max( eval(['raster_y_' num2str(iContext)]) );
        hold off
    end
    
end %ihalf


% Finish plot settings

linkaxes(hs,'x')
xlim([-Durations(1) Durations(2)])

subplot(hs(1));
set(gca,'xtick',[],'ytick',[])

subplot(hs(2))
set(gca,'ylim',[0 add_y+1],'xtick',[])

subplot(hs(3)); hold on
ylim([0 ymaxval])
xlabel('Time from transition (ms)')
ylabel('Spikes/s')
legend(ip,legstr)
hold off


% Add title
suptitle(TitleString)

figure(hfr); hold off


end