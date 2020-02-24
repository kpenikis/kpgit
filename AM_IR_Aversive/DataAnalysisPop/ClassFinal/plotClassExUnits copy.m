
% run MC_eachCell to line 140
dbstop in MC_eachCell.m at 140
MC_eachCell

if ~exist('CReach','var')
    fprintf('no Results table in workspace, so loading it...')
    tablesavename = sprintf('CR_%s.mat',whichCells);
    q=load(fullfile(figsavedir,tablesavename));
    CReach = q.CR;
end
if ~exist('CReach','var')
    CReach = CR;
end
fprintf(' done.\n')

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%% ~~~~ d'/PC vals for each stimulus

% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);


[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');

% UnitData(theseCells(iRS(iSUdps(6))))

thisUn = 56;
unitstring = sprintf('%s_%i',UnitData(theseCells(iRS(iSUdps(thisUn)))).Session,UnitData(theseCells(iRS(iSUdps(thisUn)))).Clu);


% Plot rasters
hfr=figure; 
suptitle(unitstring)
set(hfr,'Position',widescreen)

for ist = 1:size(CTTS,4)
    
    subplot(1,size(CTTS,4),ist)
    hold on
    
    for it=1:20
        sp = find(Cell_Time_Trial_Stim(theseCells(iRS(iSUdps(thisUn))),501:1000,it,ist));
        plot([sp; sp],it*ones(size([sp; sp]))+[-0.5; 0.5],'-','Color',[1 1 1]*0.2,'LineWidth',2)
        set(gca,'Color','none','ytick',[])
    end
    ylim([0 it+1.5])
end
print_eps_kp(hfr,fullfile(rootdir,'ExUnits',[unitstring '_rasters']))


% Plot PSTH
hfp=figure; 
suptitle(unitstring)
set(hfp,'Position',widescreen)

for ist = 1:size(CTTS,4)
    
    subplot(4,size(CTTS,4),ist)
    
    plot(mean(CTTS(iRS(iSUdps(thisUn)),:,:,ist),3,'omitnan')*1000,'-','Color',[1 1 1]*0.2,'LineWidth',2)
    set(gca,'Color','none')
    ylim([0 50])
end

print_eps_kp(hfp,fullfile(rootdir,'ExUnits',[unitstring '_psths']))


% Conf Mat

ConfMat = mean(CReach(iRS(iSUdps(thisUn)),:).Results{:},3);
muPC    = mean(diag(ConfMat))*100;
dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);

ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);

% Plot
hfc = figure;
imagesc(ConfMat)
axis square
caxis([-1 1])
%         cmocean('ice') %ice
cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
set(gca,'xtick',[],'ytick',[])
box off
axis off
title(sprintf('%s: %0.1f%%, d''=%0.2f',unitstring,muPC,dprime))

print_eps_kp(hfc,fullfile(rootdir,'ExUnits',[unitstring '_confmat']))




keyboard


