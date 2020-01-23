
% Taking away FR information decreases Context classification (each: small but sig; avg: ?) 

q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/Full/each/CR_each.mat');
CR_C_Full = q.CR;
q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/FullNorm/each/CR_each.mat');
CR_C_Temp = q.CR;

% Load Unit data files
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);


%% With and without FR information

figure;
plot(CR_C_Full.dprime,CR_C_Temp.dprime,'k.')
hold on
plot([0 3],[0 3])
axis square
title('each segment')

signrank(CR_C_Full.dprime, CR_C_Temp.dprime)
median(CR_C_Full.dprime - CR_C_Temp.dprime)


% Compute averages for each cell
filtCells = unique(CR_C_Full.iC);

dpHat_full    = nan(numel(filtCells),1);
dpHat_full_RS = nan(numel(filtCells),1);
dpHat_full_NS = nan(numel(filtCells),1);
for iiu = 1:numel(filtCells)
    icr = CR_C_Full.iC==filtCells(iiu);
    dpHat_full(iiu) = mean(CR_C_Full.dprime(icr));
    if ismember(filtCells(iiu),iRS)
        dpHat_full_RS(iiu) = mean(CR_C_Full.dprime(icr));
    elseif ismember(filtCells(iiu),iNS)
        dpHat_full_NS(iiu) = mean(CR_C_Full.dprime(icr));
    end
end

% % Distribution of d' RS / NS
% % No difference
% figure;
% plotSpread(dpHat_full_RS,'xValues',1,'showMM',3)
% hold on
% plotSpread(dpHat_full_NS,'xValues',2,'showMM',3)


filtCells = unique(CR_C_Temp.iC);

dpHat_temp    = nan(numel(filtCells),1);
dpHat_temp_RS = nan(numel(filtCells),1);
dpHat_temp_NS = nan(numel(filtCells),1);
for iiu = 1:numel(filtCells)
    icr = CR_C_Temp.iC==filtCells(iiu);
    dpHat_temp(iiu) = mean(CR_C_Temp.dprime(icr));
    if ismember(filtCells(iiu),iRS)
        dpHat_temp_RS(iiu) = mean(CR_C_Temp.dprime(icr));
    elseif ismember(filtCells(iiu),iNS)
        dpHat_temp_NS(iiu) = mean(CR_C_Temp.dprime(icr));
    end
end

% % Distribution of d' RS / NS
% % No difference
% figure;
% plotSpread(dpHat_temp_RS,'xValues',1,'showMM',3)
% hold on
% plotSpread(dpHat_temp_NS,'xValues',2,'showMM',3)


% PAIRWISE

% All cells
figure;
plot(dpHat_full,dpHat_temp,'k.')
hold on
plot([0 3],[0 3])
axis square
title('Avg')

signrank(dpHat_full, dpHat_temp)
median(dpHat_full - dpHat_temp)

% RS cells 
figure;
plot([0 2],[0 2])
hold on
plot(dpHat_full_RS,dpHat_temp_RS,'k.')
plot(dpHat_full_NS,dpHat_temp_NS,'r.')
axis square
title('Avg')

signrank(dpHat_full_RS, dpHat_temp_RS)
median(dpHat_full_RS - dpHat_temp_RS,'omitnan')



%% 

dpContext = sort(dpHat_temp,'descend');

figure;
plot(1:numel(dpContext),dpContext,'.k')
hold on

q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/Speech/Full/each/CR_each.mat');
CReShape = q.CR;

dpShape = sort(CReShape.dprime,'descend');

plot(1:numel(dpShape),dpShape,'.b')


skSh = skewness(dpShape,0)
skCn = skewness(dpContext,0)

% RS cells 
% Match units 









