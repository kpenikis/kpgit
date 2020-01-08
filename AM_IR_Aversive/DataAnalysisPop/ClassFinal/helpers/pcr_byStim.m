
% run MC_subpop to line 155
dbstop in MC_subpop.m at 155
MC_subpop


if ~exist('CR','var')
    fprintf('no Results table in workspace, so loading it..\n')
    tablesavename = sprintf('CR_v%s_%s.mat',varPar,whichCells);
    q=load(fullfile(figsavedir,tablesavename));
    CR = q.CR;
end


% Get results by stimulus
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


% Set stimulus order in terms of "uniqueness" 
switch whichStim
    case 'AC'
        StimOrder = [8 3 2 5 6 7 4 1];
    case 'Speech'
        StimOrder = [4 3 6 1 2 5 7 8];
end
newcolors = cmocean('thermal',8);


% Make figure
hf=figure; 
set(hf,'Position',widesmall)

subplot(1,2,1)
set(gca,'ColorOrder',newcolors);
hold on
plot(repmat([CR.iC]',[size(pcStim,1) 1])',pcStim(StimOrder,:)','LineWidth',3)
xlabel('N cells')
ylabel('PC')
grid on
set(gca,'Color','none')
xlim([0 max(CR.iC)])
ylim([0 1])
% legend(cellstr(num2str(StimOrder'))','Location','eastoutside')
% legend(cellstr(num2str(StimOrder'))','Location','southeast')
title([whichStim ', ' whichCells])

subplot(1,2,2)
set(gca,'ColorOrder',newcolors);
hold on
plot(repmat([CR.iC]',[size(dpStim,1) 1])',dpStim(StimOrder,:)','LineWidth',3)
xlabel('N cells')
ylabel('d''')
grid on
set(gca,'Color','none')
xlim([0 max(CR.iC)])
ylim([0 5.5])
legend(cellstr(num2str(StimOrder'))','Location','south')


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

% title([whichStim ', ' whichCells])


keyboard
print_eps_kp(hf,fullfile(figsavedir,'PC_byStim'))


