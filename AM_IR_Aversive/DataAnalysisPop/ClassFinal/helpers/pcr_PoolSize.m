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

% Calculate basic reults measures
pcStim  = nan(8,size(CR,1));
dpStim  = nan(8,size(CR,1));
Spsns   = nan(1,size(CR,1));
for inc = 1:size(CR,1)
    ConfMat = mean(CR(inc,:).Results{:},3);
    pcStim(:,inc)  = diag(ConfMat);
    dpStim(:,inc) = dp_from_ConfMat(ConfMat,0.001);
    Spsns(inc) = calculateSparseness(dpStim(:,inc));
end

% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);

% Sort units by: 
[pkFRsort,ipkFR] = rankPeakFR(Cell_Time_Trial_Stim(theseCells(iRS),AnWin,:,theseStim));
[dps,iSUdps]     = sort(CReach(iRS,:).dprime,'descend');

plotDPs = CReach(iRS(ipkFR),:).dprime;



% plot
hf=figure; 
set(hf,'Position',tallsmall)
hold on

% Get parameter values
PoolParams = unique([CR.iC CR.nC],'rows');
[nOccs,FirstCells] = histcounts(PoolParams(:,1),unique(PoolParams(:,1)));
FirstCells = FirstCells(nOccs>3)';

for FC = FirstCells
    
    % Sort and get indices
    thisCRidx = find(CR.iC==FC);
    [PoolSizes,sortCRidx] = sort(CR.nC(thisCRidx));
    sortCRidx = thisCRidx(sortCRidx);
    
    subplot(2,1,1)
    hold on
    plot(PoolSizes,CR.dprime(sortCRidx),'LineWidth',3)
    
    subplot(2,1,2)
    hold on
    plot(PoolSizes,Spsns(sortCRidx),'LineWidth',3)
    
end

ylabel('Sparseness')
xlabel('N cells in pool')
set(gca,'Color','none')


subplot(2,1,1)
ylabel('d''')
set(gca,'Color','none')


keyboard
print_eps_kp(gcf,fullfile(figsavedir,'PC_PoolSize'))


% translate N cells to N spikes





