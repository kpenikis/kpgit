function NunitsPerSess(SessType)

if nargin<1
    SessType = 'AM';
end

fn = set_paths_directories;

switch SessType
    case 'AM'
        q = load(fullfile(fn.processed,'Units'));
    case 'VS'
        q = load(fullfile(fn.processed,'UnitsVS'));
end
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Index units by session
[Sessions,~,idxSess] = unique(UnitInfo(:,1:2),'stable');

% i64 = find(strncmp(Sessions{idxSess,1},'A',1));
% i16 = find(strncmp(Sessions{idxSess,1},'W',1));
% 
% y64 = histcounts(idxSess(i64));
% y16 = histcounts(idxSess(i16));

y = histcounts(idxSess);

[Nuns,iss] = sort(y,'descend');

figure; 
subplot(1,6,1:5)
plot(Nuns,'.k','MarkerSize',20)
xlim([0 numel(Nuns)+1])
xlabel('Session')
ylabel('Number of SUs')
grid on
ylim([0 30])


subplot(1,6,6)
plot([0 1 1 0 0],[0 0 1 1 0],'b')
text(0.1,0.5,Sessions(iss(1:min(20,numel(iss))),:).Session)
set(gca,'xtick',[],'ytick',[])


print_eps_kp(gcf,fullfile(fn.figs,['Nunits_per_sess_' SessType]))
