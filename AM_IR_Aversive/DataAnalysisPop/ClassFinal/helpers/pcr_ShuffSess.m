
CellMin = 2;

CRplot = CR(CR.iC>CellMin,:);

% Shuffled vs sim trials
figure;
plot([0 3],[0 3],'-k')
hold on
scatter([CRplot.dprime(1:2:end)],[CRplot.dprime(2:2:end)],20.*CRplot.iC(1:2:end),'k','filled')
xlabel('d'' simultaneous trials')
ylabel('d'' shuffled trials')
axis square

print_eps_kp(gcf,fullfile(figsavedir,mastertablesavename))

