

whichStim  = 'Speech';
whichCells = 'dpRank_RS';

% Data settings
fn = set_paths_directories('','',1);

switch whichStim
    case 'AC'    
        keyboard
        
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % SHAPE: Load SU results
        datadir = fullfile(fn.figs,'ClassAM','AC','Full','each');
        q=load(fullfile(datadir,'CR_each.mat'));
        eCR_F = q.CR;
        clear q
        
        
        
        % Also get CTTS
        [~,theseCells] = recallDataParams('AC','each');
        
        savedir = fullfile(fn.figs,'ClassContext','AC','CompareShape');
        
    case 'Speech'        
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % SHAPE: Load SU results
        
        % Full
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','Full','each','CR_each.mat'));
        eCR_F = q.CR;
        clear q
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','Full',whichCells,['CR_vFull_' whichCells '.mat']));
        CR_F = q.CR;  
        clear q
        
        % Temp
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyTemp','each','CR_each.mat'));
        eCR_T = q.CR;
        clear q
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyTemp',whichCells,['CR_vOnlyTemp_' whichCells '.mat']));
        CR_T = q.CR;
        clear q
        
        % Rate
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyRate','each','CR_each.mat'));
        eCR_R = q.CR;
        clear q
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyRate',whichCells,['CR_vOnlyRate_' whichCells '.mat']));
        CR_R = q.CR;
        clear q
        
        
        % Also get CTTS
        [CTTS,theseCells] = recallDataParams('Speech','each');
        
        
        % Load Response Duration estimation
        load(fullfile(fn.figs,'ClassContext',whichStim,'RawData','avgPropUp.mat'));        
        
        savedir = fullfile(fn.figs,'ClassResults',whichStim,'Pooling');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widescreen = [1 scrsz(4)/3 scrsz(3) scrsz(4)/3];




%% For each Class Type, compare Shape and Context SU distributions

hf=figure; 
set(hf,'Position',widescreen)


% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);

% Sort units by: 
[pkFRsort,ipkFR] = rankPeakFR(CTTS(iRS,:,:,:));


% ~~~~~~~ Full ~~~~~~~

[dps,iSUdps]     = sort(eCR_F(iRS,:).dprime,'descend');

switch whichCells
    case 'dpRank_RS'
        plotDPs = eCR_F(iRS(iSUdps),:).dprime;
end

subplot(1,3,1)
plot(1:length(plotDPs),plotDPs,'.k','MarkerSize',15)
hold on
plot(CR_F.iC,CR_F.dprime,'.m','MarkerSize',25)
xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(plotDPs),6)))
xlim([0 numel(iRS)])
ylim([-0.1 4])
title('Full')


% ~~~~~~~ Temp ~~~~~~~

[dps,iSUdps]     = sort(eCR_T(iRS,:).dprime,'descend');

switch whichCells
    case 'dpRank_RS'
        plotDPs = eCR_T(iRS(iSUdps),:).dprime;
end

subplot(1,3,2)
plot(1:length(plotDPs),plotDPs,'.k','MarkerSize',15)
hold on
plot(CR_T.iC,CR_T.dprime,'.m','MarkerSize',25)
xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(plotDPs),6)))
xlim([0 numel(iRS)])
ylim([-0.1 4])
title('Temp')


% ~~~~~~~ Rate ~~~~~~~

[dps,iSUdps]     = sort(eCR_R(iRS,:).dprime,'descend');

switch whichCells
    case 'dpRank_RS'
        plotDPs = eCR_R(iRS(iSUdps),:).dprime;
end

subplot(1,3,3)
plot(1:length(plotDPs),plotDPs,'.k','MarkerSize',15)
hold on
plot(CR_R.iC,CR_R.dprime,'.m','MarkerSize',25)
xlabel('Cell N')
ylabel('d'' across stimuli')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(plotDPs),6)))
xlim([0 numel(iRS)])
ylim([-0.1 4])
title('Rate')



keyboard
print_eps_kp(gcf,fullfile(savedir,'SUvPools_15cell'))






%%  avgPropUp(iRS)


hf=figure; 
set(hf,'Position',widescreen)


% ~~~~~~~ Full ~~~~~~~
subplot(1,3,1)
plot(avgPropUp(iRS),eCR_F(iRS,:).dprime,'.k','MarkerSize',15)

xlabel('Duration proportion')
ylabel('d''')
box off
set(gca,'Color','none')
xlim([0 1])
ylim([-0.1 3])
title('Full')


% ~~~~~~~ Temp ~~~~~~~
subplot(1,3,2)
plot(avgPropUp(iRS),eCR_T(iRS,:).dprime,'.k','MarkerSize',15)

xlabel('Duration proportion')
ylabel('d''')
box off
set(gca,'Color','none')
xlim([0 1])
ylim([-0.1 3])
title('Temp')


% ~~~~~~~ Rate ~~~~~~~
subplot(1,3,3)
plot(avgPropUp(iRS),eCR_R(iRS,:).dprime,'.k','MarkerSize',15)

xlabel('Duration proportion')
ylabel('d''')
box off
set(gca,'Color','none')
xlim([0 1])
ylim([-0.1 1])
title('Rate')


keyboard
print_eps_kp(gcf,fullfile(savedir,'SUvAvgPropUp'))




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

plotDPs = CReach(iRS(iSUdps),:).dprime;


% plot
hf=figure; 
set(hf,'Position',widescreen)

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


