if ~exist('CR','var')
    % run MC_subpop to line 159
    dbstop in MC_subpop.m at 159
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


%% d', cumulative N spikes

iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
[dps,iSUdps] = sort(CReach(iRS,:).dprime,'ascend');

% idp1 = find(dps>1,1,'first');
% pBest_dps = sum(dps>1)/numel(dps);

idp1 = find(dps<prctile(dps,94),1,'last');


avgNspks = mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2),4);
pBest_nSpk = sum(avgNspks(iSUdps(idp1:end)))/sum(avgNspks);
pRest_nSpk = 1-pBest_nSpk;


ymax = 2.5;
figure; 

yyaxis right
plot(1:numel(dps),dps,'.')
hold on
plot([idp1 idp1],[-0.1 ymax],'-','Color',0.7*[1 1 1])
ylabel('dprime')
xlim([0 numel(dps)+1])
ylim([-0.1 ymax])

yyaxis left
plot(1:numel(dps),cumsum(avgNspks(iSUdps))./sum(avgNspks(iSUdps)),'-','LineWidth',3)
hold on
plot([0 idp1],[pRest_nSpk pRest_nSpk],'-','Color',0.7*[1 1 1])
ylabel('cumulative proportion of spikes')
ylim([0 1])
set(gca,'Color','none','ytick',0:0.1:1,'xtick',round((0.2:0.2:0.8)*numel(dps)))

xlim([0 numel(dps)+1])
xlabel('Ranked cells')
title(sprintf('%0.1f of spikes come from non-Best cells',pRest_nSpk*100))

print_eps_kp(gcf,fullfile(rootdir,whichStim,'SU_dp_propSpikes'))



%% compare d' to VS 

SigVS = find(~cellfun(@isempty,{UnitData(theseCells).iBMF_VS}));
NSVS  = find(cellfun(@isempty,{UnitData(theseCells).iBMF_VS}));

CReach.dprime(NSVS)


dpSg=[]; dpNS=[];
hf1=figure; 
hf2=figure; 
for icr = 1:numel(theseCells)
    
    VSdata = UnitData(theseCells(icr)).VSdata_spk(:,1:6);
    
    for ist = 2:6
        figure(hf1); hold on
        if VSdata(2,ist)>13.1
            plot(VSdata(1,ist),CReach.dprime(icr),'.b','MarkerSize',15)
            dpSg = [dpSg; CReach.dprime(icr)];
        else
            plot(VSdata(1,ist),CReach.dprime(icr),'ok','MarkerSize',3)
            dpNS = [dpNS; CReach.dprime(icr)];
        end
        figure(hf2); 
        subplot(2,3,ist); hold on
        plot(VSdata(2,ist),CReach.dprime(icr),'.k','MarkerSize',5)
    end
end
figure(hf1); hold on
xlabel('Vector strength')
ylabel('d''')
figure(hf2); hold on
xlabel('Rayleigh Statistic')
ylabel('d''')


figure;
histogram(dpNS,-0.2:0.1:2.6)
hold on
histogram(dpSg,-0.2:0.1:2.6)


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


