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
    q = load(fullfile(rootdir,whichStim,varPar,'each','CR_each.mat'));
    CReach = q.CR;
    clear q
end



pcStim  = nan(8,size(CR,1));
dpStim  = nan(8,size(CR,1));
for inc = 1:size(CR,1)
    
    ConfMat = mean(CR(inc,:).Results{:},3);
    
    pcStim(:,inc)  = diag(ConfMat);
    
    ConfMat(ConfMat==0) = 0.01;
    ConfMat(ConfMat==1) = 0.99;
    for ist = 1:size(ConfMat,1)
        othSt = 1:size(ConfMat,1);
        dpStim(ist,inc) =  norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(othSt~=ist,ist)),0,1);
    end
end


% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');


hf=figure; 
set(hf,'Position',widesmall)

subplot(1,2,1)
hold on
plot([[CR.iC]'; [CR.iC]'],[min(pcStim,[],1); max(pcStim,[],1)],'-r','LineWidth',2)
plot([CR.iC],median(pcStim,1),'.r','MarkerSize',20)

xlabel('N cells')
ylabel('min, max, median PC across stimuli')
grid on
set(gca,'Color','none')
xlim([0 numel(iRS)])
ylim([0 1])
title([whichStim ', ' whichCells])


subplot(1,2,2)
plot(1:length(dps),dps,'.k','MarkerSize',15)
hold on
% plot([[CR.iC]'; [CR.iC]'],[min(dpStim,[],1); max(dpStim,[],1)],'-r','LineWidth',2)
plot([CR.iC],median(dpStim,1),'.r','MarkerSize',20)

xlabel('N cells')
ylabel('d'' across stimuli')
grid on
set(gca,'Color','none')
xlim([0 numel(iRS)])
ylim([-0.1 5.5])


keyboard
print_eps_kp(gcf,fullfile(figsavedir,'PC_SUvPools'))


