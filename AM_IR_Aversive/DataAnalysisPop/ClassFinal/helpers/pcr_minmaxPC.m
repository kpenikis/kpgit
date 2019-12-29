
whichCells = 'Best';

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'ClassAM','AC');


% Results: adding Loudest cells
q=load(fullfile(savedir,'Full',[whichCells 'RS'],['CR_vFull_' whichCells 'RS.mat']));
CR = q.CR;

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
plot([[CR.iC]'; [CR.iC]'],[min(pcStim,[],1); max(pcStim,[],1)],'-b')
plot([CR.iC],median(pcStim,1),'.b')

xlabel('Number of cells')
ylabel('min to max PC, across stimuli')
grid on
set(gca,'Color','none')


% Add point for All cells
q=load(fullfile(savedir,'AnDur','CR_vAnDur_RS.mat'));
CRallRS = q.CR(q.CR.WinEnd==1000,:);

ConfMat = mean(CRallRS.Results{:},3);
pcStimAll = diag(ConfMat);

plot([70; 70],[min(pcStimAll,[],1); max(pcStimAll,[],1)],'-k')
plot(70,median(pcStimAll,1),'.k')

set(gca,'xtick',[[CR.iC]' 70],'xticklabel',[cellstr(num2str([CR.iC]))' 'all'])

title('Adding cells, from best SU d''')



print_eps_kp(gcf,fullfile(savedir,'Full',[whichCells 'RS'],'PC_MinMax'))




