
% run MC_subpop to line 155
dbstop in MC_subpop.m at 155
MC_subpop







if ~exist('CR','var')
    fprintf('no Results table in workspace, so loading it..\n')
    tablesavename = sprintf('CR_v%s_%s.mat',varPar,whichCells);
    q=load(fullfile(figsavedir,tablesavename));
    CR = q.CR;
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
xlim([0 max(CR.iC)])
ylim([0 1])
title([whichStim ', ' whichCells])


subplot(1,2,2)
hold on
plot([[CR.iC]'; [CR.iC]'],[min(dpStim,[],1); max(dpStim,[],1)],'-r','LineWidth',2)
plot([CR.iC],median(dpStim,1),'.r','MarkerSize',20)

xlabel('N cells')
ylabel('min, max, median d'' across stimuli')
grid on
set(gca,'Color','none')
xlim([0 max(CR.iC)])
ylim([0 5.5])

% Add point for All cells
% q=load(fullfile(fn.figs,'ClassAM','AC','AnDur','CR_vAnDur_RS.mat'));
% CRallRS = q.CR(q.CR.WinEnd==1000,:);
% 
% ConfMat = mean(CRallRS.Results{:},3);
% pcStimAll = diag(ConfMat);
% 
% plot([70; 70],[min(pcStimAll,[],1); max(pcStimAll,[],1)],'-k','LineWidth',4)
% plot(70,median(pcStimAll,1),'.k','MarkerSize',40)


% set(gca,'xtick',[CR.iC]','xticklabel',cellstr(num2str([CR.iC]))')

% suptitle([whichStim ', ' whichCells])


keyboard
print_eps_kp(gcf,fullfile(figsavedir,'PC_MinMax'))


