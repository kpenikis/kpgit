function subplot_mtx(stim,data,METRIC,lL)
% subplot_mtx(stim,data,METRIC)
%   Plots data in current axes as a heatmap, snapped to specified
%   categories for ease of viewing.
%

% Set color categories depending on metric
C=nan(size(data));
switch METRIC
    case 'FF-Pds'
        b1 = 0.4;
        b2 = 0.8;
        b3 = 1.2;
        b4 = 1.6;
    case 'VS-Pds'
        b1 = 0.15;
        b2 = 0.3;
        b3 = 0.45;
        b4 = 0.6;
    case 'RS-Pds'
        b1 = 6.9;
        b2 = 13.8;
        b3 = 20.7;
        b4 = 34.5;
end
C(data<b1)             = 1;
C(data>=b1 & data<b2)  = 2;
C(data>=b2 & data<b3)  = 3;
C(data>=b3 & data<b4)  = 4;
C(data>=b4)            = 5;

% Create color map
% map = [ 0.4  0.8  1     ;...
%         0.48 0.7  0.8   ;...
%         0.58 0.58 0.58  ;...
%         0.63 0.5  0.4   ;...
%         0.7  0.4 0.2   ];
map = [ 0.7  0.4  0.2  ;...
        0.85 0.7  0.6  ;...
        0.9  0.9  0.9  ;...
        0.55 0.55 0.55 ;...
        0.2  0.2  0.2  ];
colormap(flipud(map))

% Plot data
image(C)


% Create labels
if stim(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(stim));
else
    jitter_labels = [stim.jitter];
end
ylab = cell(1,numel(stim));
for ii = 1:numel(stim)
    
    ylab(ii) = deal({jitter_labels{ii}});
    
    if strcmp({stim(ii).behaving},'D')
        ylab{ii} = sprintf('%s: D',ylab{ii});
    end
    if strcmp({stim(ii).behaving},'A')
        ylab{ii} = sprintf('%s: A',ylab{ii});
    end
end

% Set labels
set(gca,'YTick',1:numel(stim),'YTickLabel',ylab,'TickLabelInterpreter','none')

% Other figure properties
title([num2str(stim(1).AMdepth*100) '% depth'])
xlim([0.5 size(data,2)+0.5])
ylim([0.5 size(data,1)+0.5])

% Color bar
h=colorbar; 
h.Ticks = [1 2 3 4 5];
h.TickLabels = {'0' num2str(b1) num2str(b2) num2str(b3) num2str(b4)};
ylabel(h, METRIC,'FontSize',14)

% Plot horizontal lines between conditions
ylms = get(gca,'YLim');
xlms = get(gca,'XLim');
plot(repmat(xlms,diff(ylms)-1,1)',[(ylms(1)+1):(ylms(2)-1); (ylms(1)+1):(ylms(2)-1)],'k')

% Labels on x axis only if bottom subplot
if lL
    xlabel('AM period (chronological)')
else
    set(gca,'XTickLabel',[])
end

end % subplot_mtx
