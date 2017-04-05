function subplot_bar(stim,data_mean,data_std,data_trs,baselineFR,METRIC)
% subplot_bar(stim,data_mean,data_std,data_trs,baselineFR,METRIC)
%   Plots a barplot in the current axes of the data provided.
% 


% Set labels
if stim(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(stim));
else
    jitter_labels = [stim.jitter];
end
jitters = str2double(strtok([stim.jitter],'_'));

% Get behaving datapoints
drinking = strcmp({stim.behaving},'D');
behaving = strcmp({stim.behaving},'A');

% Recalculate baseline for different conditions
if sum([drinking,behaving])>0
    baseFR = nan(size(stim));
    
    st_d = stim(drinking   &  ~behaving );
    st_p = stim( ~drinking &  ~behaving );
    st_a = stim( ~drinking &  behaving  );
    
    baseFR(drinking)              = calc_baselineFR(st_d);
    baseFR(behaving)              = calc_baselineFR(st_a);
    baseFR(~behaving & ~drinking) = calc_baselineFR(st_p);
    
else
    baseFR = baselineFR*ones(size(stim));
end

% Draw poisson line for FF
if any(strcmp(METRIC,{'VS' 'FF' 'FF-avPds'}))
    plot([0 numel(stim)+1],[1 1],'-','Color',[0.5 0.5 0.5])
elseif strcmp(METRIC,'RS')
    plot([0 numel(stim)+1],[13.8 13.8],'-','Color',[0.5 0.5 0.5])
end

% Set up vectors for creating bar plots
xbarvec = nan( 1,numel(stim));
ybarvec = nan( 1,numel(stim));
ebarvec = nan( 2,numel(stim));
xbarlab = cell(1,numel(stim));
barcols = copper(numel(stim)); %numel(unique(jitters)) because of drinking

ic=0;
for ii = 1:numel(stim)
    ic=ic+1;
    
    % Get data to plot
    xbarvec(ii) = jitters(ii);
    ybarvec(ii) = data_mean(ii);
    ebarvec(:,ii) = data_std(ii,:)';
    xbarlab(ii) = deal({jitter_labels{ii}});
    
    % Plot some things only if this is a FR plot
    if strcmp(METRIC,'FR')
        b=fill([-0.4 -0.4 0.4 0.4]+ii, [baseFR(ii) ybarvec(ii) ybarvec(ii) baseFR(ii)], barcols(ic,:));
    else
        b=bar(ii, ybarvec(ii), 'FaceColor', barcols(ic,:),'EdgeColor',barcols(ic,:));
    end
    
    switch METRIC
        % Plot FR from individual trials
        case 'FR'
            trFRs = data_trs(ii).tr;
            plot(repmat(ii,size(trFRs)), trFRs, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        case 'FF-avPds'
%             plot(repmat(ii,size(data_trs,2)), data_trs(ii,:), 'o', 'Color', barcols(ic,:),'MarkerSize',8)
    end
    
    % Special formatting by condition
    b.EdgeColor = barcols(ic,:);
    b.LineWidth = 1;
    if drinking(ii)==1
        b.FaceAlpha = 0.7;
    elseif behaving(ii)==1
        b.FaceAlpha = 0.4;
        b.EdgeColor = 'r';
    end
    if jitters(ii)==0
        b.EdgeColor = 'blue';
    end
    if stim(1).AMdepth==0
        b.EdgeColor = [0 0.4 0];
        b.LineWidth = 2;
    end
    
    % Add error bars
    errorbar(ii, ybarvec(ii), ebarvec(1,ii), ebarvec(2,ii),...
        'LineStyle','none', 'Color',[0 0 0], 'LineWidth',1 )

end

% Plot mean baseline FR and jitter=0 FR
if strcmp(METRIC,'FR')
    mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');
    plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
    plot([0 length(xbarvec)+1], [mean(baseFR) mean(baseFR)],'Color',[0.4 0.4 0.4],'LineWidth',0.2)
end

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(stim(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    switch METRIC
        case 'FR'
            ylabel('avg FR during AM')
        case 'FF'
            ylabel('avg FF during AM')
        case 'FF-avPds'
            ylabel('avg FF during AM, by periods')
        case 'VS'
            ylabel('avg Vector Strength')
        case 'RS'
            ylabel('avg Rayleigh Statistic')
        case 'standardFR'
            ylabel('avg FR during standard period')
    end
else
    set(gca,'YTickLabel',[])
end


end % subplot_bar