
% Taking away FR information decreases Context classification (each: small but sig; avg: ?) 

q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/Full/each/CR_each.mat');
CR_C_Full = q.CR;
q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/FullNorm/each/CR_each.mat');
CR_C_Temp = q.CR;

savedir = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech';

% Load Unit data files
fn = set_paths_directories('','',1);
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
UnIdx_Context = unique(CR_C_Full.iC);

dpHat_full    = nan(numel(UnIdx_Context),1);
dpHat_full_RS = nan(numel(UnIdx_Context),1);
dpHat_full_NS = nan(numel(UnIdx_Context),1);
for iiu = 1:numel(UnIdx_Context)
    icr = CR_C_Full.iC==UnIdx_Context(iiu);
    dpHat_full(iiu) = mean(CR_C_Full.dprime(icr));
    if ismember(UnIdx_Context(iiu),iRS)
        dpHat_full_RS(iiu) = mean(CR_C_Full.dprime(icr));
    elseif ismember(UnIdx_Context(iiu),iNS)
        dpHat_full_NS(iiu) = mean(CR_C_Full.dprime(icr));
    end
end

% % Distribution of d' RS / NS
% % No difference
% figure;
% plotSpread(dpHat_full_RS,'xValues',1,'showMM',3)
% hold on
% plotSpread(dpHat_full_NS,'xValues',2,'showMM',3)


UnIdx_Context = unique(CR_C_Temp.iC);

dpHat_temp    = nan(numel(UnIdx_Context),1);
dpHat_temp_RS = nan(numel(UnIdx_Context),1);
dpHat_temp_NS = nan(numel(UnIdx_Context),1);
for iiu = 1:numel(UnIdx_Context)
    icr = CR_C_Temp.iC==UnIdx_Context(iiu);
    dpHat_temp(iiu) = mean(CR_C_Temp.dprime(icr));
    if ismember(UnIdx_Context(iiu),iRS)
        dpHat_temp_RS(iiu) = mean(CR_C_Temp.dprime(icr));
    elseif ismember(UnIdx_Context(iiu),iNS)
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


% Compare within cell: d' Shape vs d' Context

% RS cells 
% Match units 

UnIdx_Context = unique(CR_C_Temp.iC);

dpHat_temp    = nan(numel(UnIdx_Context),1);
dpHat_temp_RS = nan(numel(UnIdx_Context),1);
dpHat_temp_NS = nan(numel(UnIdx_Context),1);

for iiu = 1:numel(UnIdx_Context)
    
    icr = CR_C_Temp.iC==UnIdx_Context(iiu);
    dpHat_temp(iiu) = mean(CR_C_Temp.dprime(icr));
    
    if ismember(UnIdx_Context(iiu),iRS)
        dpHat_temp_RS(iiu) = mean(CR_C_Temp.dprime(icr));
    elseif ismember(UnIdx_Context(iiu),iNS)
        dpHat_temp_NS(iiu) = mean(CR_C_Temp.dprime(icr));
    end
end


% Temporal info only
q = load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/Speech/Full/each/CR_each.mat');
CR_Shape = q.CR;

% Get iUn indices for Shape results

% Load CTTS data
[CTTS,theseCells,nUns,Dur,nStim,TrainSize,TestSize,UnitData] = recallDataParams('Speech','each');

UnIdx_Shape   = theseCells(CR_Shape.iC);
UnIdx_Context = unique(CR_C_Temp.iC);

ContextSpecialistThr = quantile(dpHat_full,.9);
ShapeSpecialistThr   = quantile(CR_Shape.dprime,.9);
ContextOkThr         = quantile(dpHat_full,.85);
ShapeOkThr          = quantile(CR_Shape.dprime,.85);

figure;
% plot([0 3],[0 3])
hold on

dpShape = [];
dpCntxt = [];
UnitIdx = [];
for ii = 1:numel(UnIdx_Context)
    iUn = UnIdx_Context(ii);
    
    if ismember(iUn,UnIdx_Shape)
        
%         icr = CR_C_Temp.iC==UnIdx_Context(ii);
        
        if ismember(iUn,iRS)
            
            plot(CR_Shape(UnIdx_Shape==iUn,:).dprime, dpHat_full(ii),'ok')
%             plot(repmat(CR_Shape(UnIdx_Shape==iUn,:).dprime,[numel(CR_C_Temp.dprime(icr)) 1]), CR_C_Temp.dprime(icr),'ob')
            dpShape = [dpShape; CR_Shape(UnIdx_Shape==iUn,:).dprime];
            dpCntxt = [dpCntxt; dpHat_full(ii)];
            UnitIdx = [UnitIdx; iUn];
            
            % Context better
            if dpHat_full(ii)>=ContextSpecialistThr && CR_Shape(UnIdx_Shape==iUn,:).dprime<ShapeOkThr %&& dpHat_full(ii)>CR_Shape(UnIdx_Shape==iUn,:).dprime
%             if dpHat_full(ii)>(2*CR_Shape(UnIdx_Shape==iUn,:).dprime) && dpHat_full(ii)>=0.5
                plot(CR_Shape(UnIdx_Shape==iUn,:).dprime, dpHat_full(ii),'og')
            end
            % Shape better
            if dpHat_full(ii)<ContextOkThr && CR_Shape(UnIdx_Shape==iUn,:).dprime>=ShapeSpecialistThr
%             if CR_Shape(UnIdx_Shape==iUn,:).dprime>(2*dpHat_full(ii)) && CR_Shape(UnIdx_Shape==iUn,:).dprime>=0.5
                plot(CR_Shape(UnIdx_Shape==iUn,:).dprime, dpHat_full(ii),'om')
            end
            % Both good
            if dpHat_full(ii)>=ContextSpecialistThr && CR_Shape(UnIdx_Shape==iUn,:).dprime>=ShapeSpecialistThr
%             if CR_Shape(UnIdx_Shape==iUn,:).dprime>(2*dpHat_full(ii)) && CR_Shape(UnIdx_Shape==iUn,:).dprime>=0.5
                plot(CR_Shape(UnIdx_Shape==iUn,:).dprime, dpHat_full(ii),'ob')
            end
            
        elseif ismember(iUn,iNS)
%             plot(CR_Shape(UnIdx_Shape==iUn,:).dprime, dpHat_full(ii),'or')
        end
    end
end

% Find cells that are good at only one task
% ContextBetter = UnitIdx(dpCntxt>(2*dpShape) & dpCntxt>=0.5);
% ShapeBetter   = UnitIdx(dpShape>(2*dpCntxt) & dpShape>=0.5);
ContextSpec = UnitIdx(dpCntxt>=ContextSpecialistThr & dpShape<ShapeOkThr);
ShapeSpec   = UnitIdx(dpCntxt<ContextOkThr & dpShape>=ShapeSpecialistThr);
BothSpec    = UnitIdx(dpCntxt>=ContextSpecialistThr & dpShape>=ShapeSpecialistThr);


figure; hold on
for ist = 1:size(CTTS,4)
    subplot(2,4,ist); hold on
    
    for ii = 1:numel(ShapeSpec)
        plot(mean(CTTS(theseCells==ShapeSpec(ii),:,:,ist),3,'omitnan'),'m')
    end
    for ii = 1:numel(ContextSpec)
        plot(mean(CTTS(theseCells==ContextSpec(ii),:,:,ist),3,'omitnan'),'g')
    end
%     for ii = 1:numel(BothSpec)
%         plot(mean(CTTS(theseCells==BothSpec(ii),:,:,ist),3,'omitnan'),'b')
%     end
end

% Response duration estimate
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/PopulationSpeechSegments/avgPropUp.mat')
figure;
plotSpread(avgPropUp(ShapeSpec),'xValues',1,'showMM',2)
plotSpread(avgPropUp(ContextSpec),'xValues',2,'showMM',2)
plotSpread(avgPropUp(BothSpec),'xValues',3,'showMM',2)



%% POOLING 

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/PoolAll/pkFR_RS/CR_pkFR_RS.mat')
CRpoolContext=CR;

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/Speech/PoolAll/pkFR_RS/CR_vPoolAll_pkFR_RS.mat')
CRpoolShape=CR;


% d' from SU (sort by peakFR)

q=load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassContext/Speech/RawData/CTTSC_Speech_sim.mat');
CTTSC = q.Cell_Time_Trial_Stim;

% Get peak FR for each stimulus 
stim_pks=nan(numel(iRS),size(CTTSC,4),size(CTTSC,5));
for iUn = 1:numel(iRS)
    for ic = 1:size(CTTSC,5)
        stim_pks(iUn,:,ic) = permute( 1000.*max(mean(CTTSC(iRS(iUn),:,:,:,ic),3,'omitnan'),[],2) ,[3 4 2 1]);
    end
end
[pkFRsort,ipkFR] = sort(median(mean(stim_pks,3,'omitnan'),2,'omitnan'),'descend');

dpSU_full_RS = dpHat_full_RS(~isnan(dpHat_full_RS));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% d' from pools

PoolSizes = unique(CRpoolContext.nC);
xvalstarts = [1 26 65];

PoolResults = nan(numel(PoolSizes),numel(xvalstarts),numel(unique(CRpoolContext.Seg)));

for iSeg = unique(CRpoolContext.Seg)'
    
    icr = CRpoolContext.Seg==iSeg;
    CellStarts = unique(CRpoolContext(icr,:).iC)';
    
    for cs = 1:numel(CellStarts)
        PoolResults(:,cs,iSeg) = CRpoolContext(CRpoolContext(icr,:).iC==CellStarts(cs),:).dprime;    
    end
end

% Plot 
figure;
hold on
plot( dpSU_full_RS(ipkFR),'.k')
for cs = 1:numel(xvalstarts)
    plot(xvalstarts(cs)+PoolSizes,mean(PoolResults(:,cs,:),3,'omitnan'),'-','LineWidth',2)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Separate plot for each Segment
%   ** still need to filter SU data

PoolSizes = unique(CRpoolContext.nC);

for iSeg = unique(CRpoolContext.Seg)'
    
    figure;
    hold on
    plot( dpSU_full_RS(ipkFR),'.k')
    
    icr = find(CRpoolContext.Seg==iSeg);
    CellStarts = unique(CRpoolContext(icr,:).iC)';
    
    for cs = 1:numel(CellStarts)
        plot(PoolSizes+CellStarts(cs)*ones(size( CRpoolContext(CRpoolContext(icr,:).iC==CellStarts(cs),:).dprime)),...
            CRpoolContext(CRpoolContext(icr,:).iC==CellStarts(cs),:).dprime,'LineWidth',2);    
    end
    
    title(CRpoolContext(icr(1),:).SegStr{:})
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Running pool of size 15
xvalstarts = ceil(numel(dpSU_full_RS)*[0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.7]);
PoolSizes = 15; 
PoolResults = nan(numel(PoolSizes),numel(xvalstarts),numel(unique(CRpoolContext.Seg)));

for iSeg = unique(CRpoolContext.Seg)'
    
    icr = CRpoolContext.Seg==iSeg;
    CellStarts = unique(CRpoolContext(icr,:).iC)';
    
    for cs = 1:numel(CellStarts)
        
        PoolResults(:,cs,iSeg) = CRpoolContext( CRpoolContext(icr,:).iC==CellStarts(cs) & CRpoolContext(icr,:).nC==PoolSizes ,:).dprime;    
    end
end

figure;
hold on
plot( dpSU_full_RS(ipkFR),'.k')
plot( xvalstarts+PoolSizes/2, mean(PoolResults,3,'omitnan'),'LineWidth',2)



% Compare to Shape

ImprvScore_Context = [];

figure;
plot([0 4],[0 4])
hold on
for icr = 1:size(CRpoolContext,1)
    plot( max(CRpoolContext(icr,:).SUdps{:}), CRpoolContext(icr,:).dprime,'ob')
%     ImprvScore_Context = [ImprvScore_Context; (CRpoolContext(icr,:).dprime-max(CRpoolContext(icr,:).SUdps{:}))/max(CRpoolContext(icr,:).SUdps{:})];
    ImprvScore_Context = [ImprvScore_Context; (CRpoolContext(icr,:).dprime-max(CRpoolContext(icr,:).SUdps{:}))/max(CRpoolContext(icr,:).SUdps{:})];
end


ImprvScore_Shape = [];

for icr = 1:size(CRpoolShape,1)
    plot( max(CRpoolShape(icr,:).SUdps{:}), CRpoolShape(icr,:).dprime,'ok')
%     ImprvScore_Shape = [ImprvScore_Shape; (CRpoolShape(icr,:).dprime-max(CRpoolShape(icr,:).SUdps{:}))/max(CRpoolShape(icr,:).SUdps{:})];
    ImprvScore_Shape = [ImprvScore_Shape; (CRpoolShape(icr,:).dprime-max(CRpoolShape(icr,:).SUdps{:}))/max(CRpoolShape(icr,:).SUdps{:})];
end

figure;
plotSpread(ImprvScore_Shape,'xValues',1,'showMM',2)
plotSpread(ImprvScore_Context,'xValues',2,'showMM',2)




