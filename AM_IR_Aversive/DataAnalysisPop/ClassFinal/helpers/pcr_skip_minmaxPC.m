
pcStim  = nan(8,size(CR,1));
dpStim  = nan(8,size(CR,1));
for inc = 1:size(CR,1)
    
    ConfMat = mean(CR(inc,:).Results{:},3);
    
    pcStim(:,inc)  = diag(ConfMat);
    
    ConfMat(ConfMat==0) = 0.0001;
    ConfMat(ConfMat==1) = 0.9999;
    for ist = 1:size(ConfMat,1)
        othSt = 1:size(ConfMat,1);
        dpStim(ist,inc) =  norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(ist,othSt~=ist)),0,1);
    end
end


hsf1=figure; hold on
plot([[CR.iC]'; [CR.iC]'],[min(pcStim,[],1); max(pcStim,[],1)],'-b','LineWidth',4)
plot([CR.iC],median(pcStim,1),'.b','MarkerSize',40)

xlabel('N cells skipped')
ylabel('min, max, and median PC across stimuli')
grid on
set(gca,'Color','none')
ylim([0 1])

% Add point for All cells
q=load(fullfile(fn.figs,'ClassAM','AC','AnDur','CR_vAnDur_RS.mat'));
CRallRS = q.CR(q.CR.WinEnd==1000,:);

ConfMat = mean(CRallRS.Results{:},3);
pcStimAll = diag(ConfMat);

plot([70; 70],[min(pcStimAll,[],1); max(pcStimAll,[],1)],'-k','LineWidth',4)
plot(70,median(pcStimAll,1),'.k','MarkerSize',40)

set(gca,'xtick',[[CR.iC]' 70],'xticklabel',[cellstr(num2str([CR.iC]))' 'all'])

title('Removing cells, from best SU d''')


print_eps_kp(gcf,fullfile(savedir,'PC_MinMax'))



