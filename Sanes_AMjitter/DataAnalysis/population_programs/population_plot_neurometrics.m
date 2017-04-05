function population_plot_neurometrics(dp_struct,stim,plotOptions,ALLcolors,hS,xlimits)

% Get unique conditions for discriminated variable
for ic = 1:numel(dp_struct.data)
    
    % Set colors according to this condition
    this_condition = stim{ic,1};
    
    plotOptions.dataColor = ALLcolors( strcmp(strtok(this_condition,'_'),strtok(plotOptions.colSelect,'_')) , : );
    plotOptions.lineColor = plotOptions.dataColor;
    
    % Plot data
    plot_neurometric( dp_struct.data(ic), dp_struct.fitdata(ic), plotOptions, hS, xlimits, plotOptions.dataColor );
    
end

end