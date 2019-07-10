

1316
1326
1027


AMrates = [2 4 8 16 32];
xvals = 1:5;

thesecolors = cmocean('phase',5);

htc=figure;
hold on





ic = 4;
y_vals = FR_1027;
y_errs = std_1027;

% Plot baseline FR

plot([0 max(xvals)+1],[y_vals(9) y_vals(9)],'--','Color',thesecolors(ic,:),'LineWidth',3)

% Plot periodic stimuli
plot(1:5,y_vals(2:6),'Color',thesecolors(ic,:),'LineWidth',2)
for ir = 1:5
    plot([xvals(ir) xvals(ir)], y_vals(ir+1) + y_errs(ir+1)*[-1 1], ...
        '-','Color',thesecolors(ic,:),'LineWidth',4)
    plot(xvals(ir),y_vals(ir+1), 'o','MarkerSize',20,...
        'MarkerFaceColor',thesecolors(ic,:),'MarkerEdgeColor','none')
end







% Finish formatting
set(gca,'XTick',min(xvals):max(xvals),...
    'XTickLabel',AMrates,...
    'TickLabelInterpreter','none')
xtickangle(45)

xlim([-1+min(xvals) max(xvals)+1])
ylim([0 60])
ylabel('FR')
