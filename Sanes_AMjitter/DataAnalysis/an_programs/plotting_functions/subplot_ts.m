function subplot_ts(stim,r_mean,p_mean,shift)
% subplot_ts(stim,r_mean,p_mean,shift)
%   Plots a line plot, with dot at max.

% Set labels
if stim(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(stim));
else
    jitter_labels = [stim.jitter];
end
% ylab(ii) = deal({jitter_labels{ii}});

leg_str = cell(1,numel(stim));
for ii = 1:numel(stim)
    
    leg_str(ii) = deal({jitter_labels{ii}});
    
    if strcmp({stim(ii).behaving},'D')
        leg_str{ii} = sprintf('%s: D',leg_str{ii});
    end
    if strcmp({stim(ii).behaving},'A')
        leg_str{ii} = sprintf('%s: A',leg_str{ii});
    end
end

jitters = str2double(strtok([stim.jitter],'_'));

% Get behaving datapoints
drinking = strcmp({stim.behaving},'D');
behaving = strcmp({stim.behaving},'A');

barcols = copper(numel(stim)); %numel(unique(jitters)) because of drinking

plot([0 0],[-0.15 0.15],'--k')

for is = 1:numel(stim)
    [rmax,im] = max(r_mean(is,:));
    if drinking(is)==1
        hp(is)=plot(shift,r_mean(is,:),'--','Color',barcols(is,:),'LineWidth',3);
        plot(shift(im),rmax,'o','MarkerEdgeColor','b','MarkerFaceColor',barcols(is,:),'MarkerSize',15)
    elseif behaving(is)==1
        hp(is)=plot(shift,r_mean(is,:),'--','Color',barcols(is,:),'LineWidth',3);
        plot(shift(im),rmax,'o','MarkerEdgeColor','r','MarkerFaceColor',barcols(is,:),'MarkerSize',15)
    else
        hp(is)=plot(shift,r_mean(is,:),'Color',barcols(is,:),'LineWidth',3);
        plot(shift(im),rmax,'o','MarkerEdgeColor',barcols(is,:),'MarkerFaceColor',barcols(is,:),'MarkerSize',15)
    end
end
legend(hp,leg_str,'Interpreter','none','Location','southwest')

% Set axis properties
xlim([min(shift) max(shift)]) 
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
xlabel({'Shift amount (ms)'; 'positive: waveform shifts later'})
ylabel('avg Pearson R, corr between stim and individual trials')
end
title([num2str(stim(1).AMdepth*100) '% depth'])


end