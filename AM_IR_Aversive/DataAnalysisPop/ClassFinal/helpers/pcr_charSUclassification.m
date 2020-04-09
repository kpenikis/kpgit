
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


%% ~~~~ d'/PC vals for each stimulus

% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<0.43);


[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');

iRS(iSUdps)

maxdps = max(dpStim,[],2);
[maxdps_pl,imaxdps] = sort(maxdps,'descend');

figure;
plot(1:length(dps),dps,'.k','MarkerSize',15)
hold on
plot(1:numel(maxdps_pl),maxdps_pl,'.m','MarkerSize',15)
xlim([0 length(dps)])
ylim([-0.1 5])
set(gca,'Color','none')
grid on

print_eps_kp(gcf,fullfile(figsavedir,'SUdps_sort'))



% figure;
% % set(gcf,'Position',widesmall)
% plot(1:length(dps),dps,'.k','MarkerSize',15)
% xlim([0 length(dps)])
% ylim([-0.05 2.5])
% set(gca,'ytick',0:0.5:2.5,'xtick',[5 10 30 50 100 150])
% grid on
% title([whichStim ' -- SU dprimes'])
% ylabel('d prime')
% xlabel('Cell number')

% keyboard
% print_eps_kp(gcf,fullfile(savedir,'cdf_SUdps'))


% For each cell, get its PC & d' for each stimulus
pcStim = nan(length(iRS),nStim);
dpStim = nan(length(iRS),nStim);
for ii = 1:length(iRS)
        
    ConfMat = mean(CReach(iRS(ii),:).Results{:},3,'omitnan');
    
    pcStim(ii,:) = diag(ConfMat)';
    
    ConfMat(ConfMat==0) = 0.01;
    ConfMat(ConfMat==1) = 0.99;
    for ist = 1:size(ConfMat,1)
        othSt = 1:size(ConfMat,1);
        dpStim(ii,ist) =  norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(othSt~=ist,ist)),0,1);
    end
end

% Sequence a little different when d' calculated this way
% [meandps,iUdp] = sort(mean(dpStim,2),'descend');


switch whichStim
    case 'AC'
        StimOrder = [8 3 2 5 6 7 4 1];
    case 'Speech'
        StimOrder = [4 3 6 1 2 5 7 8];
end
newcolors = cmocean('thermal',8);


close all
for ist = StimOrder
    
    hf(ist)=figure;
    set(gcf,'Position',widesmall./[1 1 2 1])
    
    hold on
    
    % Add this stimulus dp for each cell
    plot([1:length(dps); 1:length(dps)],[zeros(length(iSUdps),1) dpStim(iSUdps,ist)]',...
        '-','Color',newcolors(ist==StimOrder,:),'LineWidth',2)
    
    plot(1:length(dps),dps,'.k','MarkerSize',10)
    
    xlim([0 length(dps)+1])
    ylim([-0.05 2.5])
    set(gca,'ytick',0:0.5:2.5,'xtick',[10 30 50 100 150])%,'Color','none')
    grid on
    title(['SU dprimes -- st# ' num2str(ist)])
    ylabel('d prime')
    xlabel('Cell number')
    
    print_eps_kp(hf(ist),fullfile(figsavedir,['SUdps_st' num2str(ist)]))
end




%% Compare d' to FF 



%%
% nSpks  = (mean(sum(mean(CTTS,3,'omitnan'),2,'omitnan'),4));
% nSpkRS = (mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2,'omitnan'),4));
% nSpkNS = (mean(sum(mean(CTTS(iNS,:,:,:),3,'omitnan'),2,'omitnan'),4));


figure;
subplot(1,4,1)
hold on

% Manually make boxplots
q5  = quantile(CReach.dprime,0.05);
q25 = quantile(CReach.dprime,0.25);
q75 = quantile(CReach.dprime,0.75);
q95 = quantile(CReach.dprime,0.95);

plot([1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
fill(1+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')

ylabel('dprime')
set(gca,'Color','none')
box off
ylim([-0.5 2.5])
plotSpread(CReach.dprime,'distributionIdx',ones(size(CReach.dprime)),'distributionColors','k','showMM',3)


subplot(1,4,2:4)
%         plot(nSpks,CReach.dprime,'k.')
hold on
plot(nSpkRS,CReach.dprime(iRS),'r.')
plot(nSpkNS,CReach.dprime(iNS),'b.')
ylim([-0.5 2.5])
set(gca,'ytick',[],'Color','none')
box off
xlabel('N spks')
%         title('SU dprimes vs mean N spks per stim')


% Stats
[r_RS,p_RS] = corr(nSpkRS,CReach.dprime(iRS),'Type','Spearman');
[r_NS,p_NS] = corr(nSpkNS,CReach.dprime(iNS),'Type','Spearman');
[r,p]       = corr(nSpks,CReach.dprime,'Type','Spearman');

text(40,-0.3,sprintf('Spearman r=%0.2f, p=%0.3e',r,p))


savename = sprintf('SU_dps_%s_%s',varPar,whichStim);

keyboard

print_eps_kp(gcf,fullfile(savedir,savename))

