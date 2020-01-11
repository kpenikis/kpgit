% close all
crAmt = 0.001;

whichStim = 'Speech'; %finish generalizing


% Data settings
fn = set_paths_directories('','',1);

switch whichStim
    case 'AC'
        datadir = fullfile(fn.figs,'ClassAM',whichStim,'Full','each');
        savedir = fullfile(fn.figs,'ClassAM',whichStim,'Full','each','byRespChar');
        
        q = load(fullfile(fn.processed,'Units'));
    case 'Speech'
        datadir = fullfile(fn.figs,'ClassSpeech',whichStim,'Full','each');
        savedir = fullfile(fn.figs,'ClassSpeech',whichStim,'Full','each','byRespChar');
        
        q = load(fullfile(fn.processed,'UnitsVS'));
end

% Unit data files
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load SU data
tabsavename = sprintf('CR_%s.mat','each');
q=load(fullfile(datadir,tabsavename));
CR_SU = q.CR;
clear q

% Also get CTTS
[CTTS,theseCells,nUns,Dur,nStim,TrainSize,TestSize,UnitData] = recallDataParams(whichStim,'each');


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

% keyboard

%%

% iSigUns = identifyResponsiveUnits(UnitData(theseCells));
% iNSuns  = 1:length(theseCells);
% iNSuns(iSigUns) = [];
% 
% figure;
% plotSpread(CR_SU(iSigUns,:).dprime,'distributionIdx',ones(size(CR_SU(iSigUns,:).dprime)))
% hold on
% plotSpread(CR_SU(iNSuns,:).dprime,'distributionIdx',1+ones(size(CR_SU(iNSuns,:).dprime)))


%%

% Plot peak FR by d'

[pkFRsort,ipkFR] = rankPeakFR(CTTS);

figure; hold on

yyaxis right
plot(1:length(ipkFR),CR_SU.dprime(ipkFR),'.')
ylabel('d'', median across stimuli')
ylim([-0.2 3.5])

yyaxis left
plot(1:length(ipkFR),pkFRsort,'.')
ylabel('peak FR, median across stimuli')
xlim([0 numel(pkFRsort)+1])

print_eps_kp(gcf,fullfile(datadir,'dp_peakFR'))

keyboard

% Correlation: overall, avg across stim
figure;
plot(median(stim_pks,2),CR_SU.dprime,'.')
hold on
plot(max(stim_pks,[],2),CR_SU.dprime,'.','MarkerSize',20)

xlim([5 1000])
ylim([-0.2 2.5])
ylim([0.0001 10])
set(gca,'xscale','log','Color','none')
xlabel('max (+median) of peak FR across stimuli')
ylabel('SU d''')

[rmed,pmed] = corr(median(stim_pks,2),CR_SU.dprime);
[rmax,pmax] = corr(max(stim_pks,[],2),CR_SU.dprime);
[rmax,pmax] = corr(log(max(stim_pks,[],2)),real(log(CR_SU.dprime)));


% each stimulus (must re-calc d') 
pcStim = nan(size(CR_SU,1),8);
dpStim = nan(size(CR_SU,1),8);
for iu = 1:size(CR_SU,1)
    ConfMat = mean(CR_SU(iu,:).Results{:},3);
    pcStim(iu,:)  = diag(ConfMat)';
    dpStim(iu,:) = dp_from_ConfMat(ConfMat,0.001);
end

figure;
for ist = 1:8
    subplot(4,2,ist)
    plot(stim_pks(:,ist),dpStim(:,ist),'.')
    title(ist)
    xlim([5 1000])
    ylim([-0.2 4.5])
    set(gca,'xscale','log','Color','none')
    xlabel('peak FR')
    ylabel('SU d''')
end

% keyboard


%%

% Plot Mean vs Variance, color by d'

% Each stim 
stim_mus=nan(size(CTTS,1),8);
stim_var=nan(size(CTTS,1),8);
for iUn = 1:size(CTTS,1)
    stim_mus(iUn,:) = mean( permute(sum(CTTS(iUn,:,:,:),2),[3 4 1 2]) ,1,'omitnan');
    stim_var(iUn,:) = var(  permute(sum(CTTS(iUn,:,:,:),2),[3 4 1 2]) ,1,'omitnan');
end

figure;
plot(stim_mus,stim_var,'.')
xlim([0 1400])
ylim([0 1400])


figure;
scatter(mean(stim_mus,2),mean(stim_var,2),40,CR_SU.dprime,'filled')
xlim([0 600])
ylim([0 600])
axis square

aaa=245;













