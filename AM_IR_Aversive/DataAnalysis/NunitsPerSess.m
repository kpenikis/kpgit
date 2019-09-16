
fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Index units by session
[Sessions,~,idxSess] = unique([UnitInfo.Subject UnitInfo.Session],'stable');

y = histogram(idxSess);
[Nuns,iss] = sort(y.Values,'descend');

figure; 
plot(Nuns,'.k','MarkerSize',20)
xlabel('Session')
ylabel('Number of SUs')
grid on

print_eps_kp(gcf,fullfile(fn.figs,'Nunits_per_sess'))
