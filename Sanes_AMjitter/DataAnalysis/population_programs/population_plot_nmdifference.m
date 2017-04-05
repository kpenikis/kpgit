function population_plot_nmdifference(dp_struct,stim,plotOptions,ALLcolors,hS,xlimits)

% Get unique conditions for discriminated variable

% Get dprimes for jitter=0 or depth=0 condition
for ic = 1:numel(dp_struct.data)
    
    this_condition = stim{ic,1};
    
    if strcmp(strtok(this_condition,'_'),'0')
        baseline = dp_struct.data{ic}(:,2);
    else
        continue
    end 
end

% Plot the difference for the rest of the conditions
for ic = 1:numel(dp_struct.data)
    
    this_condition = stim{ic,1};
    
    if strcmp(strtok(this_condition,'_'),'0')
        continue
    end
        
    plotOptions.dataColor = ALLcolors( strcmp(strtok(this_condition,'_'),strtok(plotOptions.colSelect,'_')) , : );
    plotOptions.lineColor = plotOptions.dataColor;
    
    % Plot nm difference data
try
    plot(dp_struct.data{ic}(:,1),dp_struct.data{ic}(:,2) - baseline,'.-',...
        'MarkerSize',plotOptions.dataSize,...
        'Color',plotOptions.dataColor);
    hold on;
    xlim(xlimits)
    set(gca,'TickDir','out')
    box off
catch
    keyboard
end
end

end
