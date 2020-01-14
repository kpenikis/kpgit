if ~exist('CR','var')
    % run MC_subpop to line 157
    dbstop in MC_subpop.m at 157
    MC_subpop
end

% Get pooled cells results
if ~exist('CR','var')
    fprintf('no Results table in workspace, so loading it..\n')
    tablesavename = sprintf('CR_v%s_%s.mat',varPar,whichCells);
    q=load(fullfile(figsavedir,tablesavename));
    CR = q.CR;
    clear q
end

% Load SU results
if ~exist('CReach','var')
    q = load(fullfile(rootdir,whichStim,'Full','each','CR_each.mat'));
    CReach = q.CR;
    clear q
end


%%

pcStim  = nan(8,size(CR,1));
dpStim  = nan(8,size(CR,1));
dpStim2  = nan(8,size(CR,1));
for inc = 1:size(CR,1)
    
    ConfMat = mean(CR(inc,:).Results{:},3);
    
    pcStim(:,inc)  = diag(ConfMat);
    
    dpStim(:,inc) = dp_from_ConfMat(ConfMat,0.001);
    
end


% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);

% Sort units by: 
[pkFRsort,ipkFR] = rankPeakFR(CTTS(iRS,:,:,:));
[dps,iSUdps]     = sort(CReach(iRS,:).dprime,'descend');

plotDPs = CReach(iRS(ipkFR),:).dprime;


% plot
hf=figure; 
set(hf,'Position',widesmall)

subplot(1,2,1)
hold on
% plot([[CR.iC]'; [CR.iC]'],[min(pcStim,[],1); max(pcStim,[],1)],'-r','LineWidth',2)
% plot([CR.iC],median(pcStim,1),'.r','MarkerSize',20)
plot(CR.iC,CR.PC./100,'.m','MarkerSize',20)

xlabel('Cell N')
ylabel('min, max, median PC across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(plotDPs),6)))
xlim([0 numel(iRS)])
ylim([0 1])
title([whichStim ', ' whichCells])


subplot(1,2,2)
plot(1:length(plotDPs),plotDPs,'.k','MarkerSize',15)
hold on
% plot([[CR.iC]'; [CR.iC]'],[min(dpStim,[],1); max(dpStim,[],1)],'-r','LineWidth',2)
% plot([CR.iC CR.iC+CR.nC-1]',[median(dpStim,1); median(dpStim,1)],'>r','MarkerSize',20)
% plot([CR.iC CR.iC+CR.nC-1]',[median(dpStim,1); median(dpStim,1)],'-r','LineWidth',2)
% plot([CR.iC CR.iC+CR.nC-1]',[CR.dprime'; CR.dprime'],'-m','LineWidth',2)
plot(CR.iC,CR.dprime,'.m','MarkerSize',25)

xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(plotDPs),6)))
xlim([0 numel(iRS)])
ylim([-0.1 4.5])


keyboard
print_eps_kp(gcf,fullfile(figsavedir,'PC_SUvPools'))


