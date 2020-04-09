
% Clustered matrix of SU d' by stimulus
% Pairwise correlation of stimulus SU d' vectors
%   subsequently added:
% Average confusion matrices
% Boxplots of SU d' by stimulus

close all

whichClass   = 'Full';
whichStim    = 'AC';


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',14)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%% Data 

fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,whichStim);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

switch whichStim
    case 'AC'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassAM',whichStim,whichClass);
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        istlab = {'Unmod' '2hz' '4hz' '8hz' '16hz' '32hz' 'IrrA' 'IrrB'};
        
    case 'Speech'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassSpeech',whichStim,whichClass);
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        istlab = {'ber' 'I can''t' 'blab' 'I''m' 'which you' 'to be chop' 'please' 'trees'};
        
end

% Load SU data 
q=load(fullfile(datadir,'each','CR_each.mat'));
CReach = q.CR;
clear q

% Load subpop results    
q=load(fullfile(datadir,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
CR_pfr = q.CR;
clear q
q=load(fullfile(datadir,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
CR_Qpfr = q.CR;
clear q
q=load(fullfile(datadir,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
CR_dp = q.CR;
clear q

% Also get CTTS
[CTTS,theseCells] = recallDataParams(whichStim,'each',12);


%%
% RS indices
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);

% Significance check
% UnSig = bootstrap4significance(CReach(iRS,:));

keyboard

%% Average confusion matrix
%  For SU, 16 Hz often confused with Warn and 32.
%  For rand pools of 10 cells, PC of 16 Hz *not* better in relative terms
%  It does look like FA are common for Warn with SU, and not with pools


% SU
foo = cellfun(@(x) mean(x,3,'omitnan'),CReach.Results,'UniformOutput',false);
MeanConfMat_SU = mean( cat(3,foo{:}), 3);

figure;
imagesc(MeanConfMat_SU-0.125)
axis square
caxis([-0.2 0.2])
cmocean('curl','pivot',0)

% Rand pools of 10
foo = cellfun(@(x) mean(x,3,'omitnan'),CR.Results,'UniformOutput',false);
MeanConfMat_Rand = mean( cat(3,foo{:}), 3);

figure;
imagesc(MeanConfMat_Rand-0.125)
axis square
caxis([-0.2 0.2])
cmocean('curl','pivot',0)


figure;
plot(diag(MeanConfMat_SU)/mean(diag(MeanConfMat_SU)),'LineWidth',2)
hold on
plot(diag(MeanConfMat_Rand)/mean(diag(MeanConfMat_Rand)),'LineWidth',2)


% But wait! Now just cells with 16 Hz VS>0.5
iC16 = cell2mat(cellfun(@(x) x(1,5)>=0.5,{UnitData(theseCells).VSdata_spk},'UniformOutput',false))';
iC16 = cell2mat(cellfun(@(x) x(2,5)>=13.1,{UnitData(theseCells).VSdata_spk},'UniformOutput',false))';
% UnitData(theseCells(iC16)).VSdata_spk
foo = cellfun(@(x) mean(x,3,'omitnan'),CReach(iC16,:).Results,'UniformOutput',false);
MeanConfMat_SU = mean( cat(3,foo{:}), 3);

figure;
imagesc(MeanConfMat_SU-0.125)
axis square
caxis([-0.2 0.2])
cmocean('curl','pivot',0)

% still doesn't suggest that 16 is misclassified as 8 or 4 more often than
% chance (as might be expected if periods were missed)



%%

dpmaxval = 3;
nStim    = size(CReach(1,:).Results{:},1);


% For each cell, get its PC & d' for each stimulus
pcStim = nan(length(iRS),nStim);
dpStim = nan(length(iRS),nStim);
Sparss = nan(length(iRS),1);
for ii = 1:length(iRS)
        
    ConfMat = mean(CReach(iRS(ii),:).Results{:},3,'omitnan');
    
    pcStim(ii,:) = diag(ConfMat)';
    dpStim(ii,:) = dp_from_ConfMat(ConfMat,0.01);
    Sparss(ii) = calculateSparseness(dpStim(ii,:)');
end


%==========================================================================
%                     boxplots of SU d' by stimulus
%==========================================================================

hfb=figure;
set(hfb,'Position',widescreen)
hold on

for ist = 1:size(dpStim,2)
    
    % Manually make boxplots
    q5  = quantile(dpStim(:,ist),0.05);
    q25 = quantile(dpStim(:,ist),0.25);
    q75 = quantile(dpStim(:,ist),0.75);
    q95 = quantile(dpStim(:,ist),0.95);
    
    plot([ist ist],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
    fill(ist+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
        
end

P = anova1(dpStim)

plotSpread(dpStim,'distributionIdx',ist*ones(size(dpStim,1)),'distributionColors','k','showMM',3)

ylabel('dprime')
set(gca,'Color','none')
box off
% ylim([-0.5 5])



%==========================================================================
%                               SU d' 
%                    cluster cells by tuning prefs
%==========================================================================

DATA = round(dpStim,1);
idplab = sprintfc('%d',1:size(DATA,1));

hf=figure; 
set(hf,'Position',widescreen)

% Calulate linkages
Y = pdist(DATA);     % cells X time
Z = linkage(Y,'ward');

leafOrder = fliplr(optimalleaforder(Z,Y));


subplot(2,1,1)     % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'Orientation','top','Labels',idplab);%,'reorder',leafOrder);%,'Labels',idplab);
set(gca,'tickdir','out')


% Clustered by tuning curve similarity
subplot(2,1,2)
imagesc( DATA(fliplr(outperm),:)' )
caxis([-0.2 dpmaxval])

DATA_sc = DATA(fliplr(outperm),:)';
sortCells = fliplr(outperm);


%%

DATA2 = DATA_sc;

% Calulate linkages
Y = pdist(DATA2);     % cells X time
Z = linkage(Y,'ward');

leafOrder = fliplr(optimalleaforder(Z,Y));


% Plot
hf2=figure; 
set(hf2,'Position',widescreen)

subplot(7,4,[1 5])     % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'Orientation','top','Labels',istlab,'reorder',leafOrder);%,'reorder',leafOrder);%,'Labels',idplab);
set(gca,'tickdir','out')
xtickangle(45)


% Clustered by tuning curve similarity
subplot(7,4,(4*3+1):(4*6))
imagesc( DATA2(fliplr(outperm),:) )
caxis([-0.2 dpmaxval])
set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(fliplr(outperm)))
hcb=colorbar;
hcb.Label.String = 'd''';
xlabel('RS cells')


print_eps_kp(hf2,fullfile(savedir,'SUdp_stim_sorted'))

sortStim = fliplr(outperm);
DATA_sorted = DATA_sc(fliplr(outperm),:);


%% Plot same matrix when cells sorted by mean d'

[meandp,imdp] = sort(CReach.dprime(iRS),'descend');

hf3=figure; 
set(hf3,'Position',widescreen)

subplot(7,4,(4*3+1):(4*6))
imagesc( DATA(imdp,:)' )
caxis([-0.2 dpmaxval])
set(gca,'ytick',1:numel(istlab),'yticklabel',istlab)
hcb=colorbar;
hcb.Label.String = 'd''';
xlabel('RS cells')

print_eps_kp(hf3,fullfile(savedir,'SUdp_sortMeanDp'))

[~,iC_mdp] = sort(imdp);
[~,iC_cls] = sort(sortCells');

[r,p] = corr(iC_mdp,iC_cls)


%%
%==========================================================================
% stim d' compared to avg d'
%==========================================================================

hf1=figure;
set(hf1,'Position',fullscreen)

for ist = 1:size(dpStim,2)
    
    subplot(2,4,ist)
    plot(dpStim(:,ist),CReach.dprime(iRS),'ok')
    
    axis square
    set(gca,'Color','none')
    xlim([-0.5 4])
    ylim([-0.5 4])
    xlabel('d'' of this stim')
    ylabel('avg d''')
    title(istlab{ist})
    
    [r,p] = corr(dpStim(:,ist),CReach.dprime(iRS));
    text(0,3.5,sprintf('r=%0.2f',r),'Color','b','FontSize',14)
end

print_eps_kp(hf1,fullfile(savedir,'SUdp_mean-stim'))


%%
%==========================================================================
% Try to segregate rate coding and phase locking cells
%==========================================================================

figure;
plot(dpStim(:,6),max(dpStim(:,2:4),[],2),'ok')
axis square
xlim([-0.2 3.5])
ylim([-0.2 3.5])
xlabel('SU d'' 32 Hz')
ylabel('SU d'' max from 2, 4, or 8 Hz')


figure;
plot(dpStim(:,6),mean(dpStim(:,2:4),2),'ok')
axis square
xlim([-0.2 3.5])
ylim([-0.2 3.5])
xlabel('SU d'' 32 Hz')
ylabel('SU mean d'' of 2,4,8 Hz')

print_eps_kp(gcf,fullfile(savedir,'SUdp_slow_vs_32'))


%%
%==========================================================================
% Pairwise correlations betweene stimulus SU d' vectors
%==========================================================================

CorrMat = ones(nStim,nStim);
for ist1=1:nStim
    for ist2=(ist1+1):nStim
        CorrMat(ist1,ist2) = corr(dpStim(:,ist1),dpStim(:,ist2),'type','Pearson');
        CorrMat(ist2,ist1) = corr(dpStim(:,ist1),dpStim(:,ist2),'type','Pearson');
    end
end

% Periodic stimuli only
figure;
imagesc(CorrMat(2:6,2:6))
axis square
cmocean('grey')
caxis([0.5 1])
colorbar
set(gca,'xtick',1:5,'xticklabel',istlab(2:6),'ytick',1:5,'yticklabel',istlab(2:6))
title('Pairwise correlations of SU d''s')
print_eps_kp(gcf,fullfile(savedir,'SUdp_xcorrPdcstim'))


% Sort stimuli by clustering algorithm
figure;
imagesc(CorrMat(sortStim,sortStim))
axis square
cmocean('grey')
caxis([0.5 1])
colorbar
set(gca,'xticklabel',istlab(sortStim),'yticklabel',istlab(sortStim))
title('Pairwise correlations of SU d''s')

print_eps_kp(gcf,fullfile(savedir,'SUdp_xcorrstim_clustered'))



%==========================================================================
%%                           SU d' per stim
%                   correlate to response properties
%==========================================================================

%%  VS x d' 

alfa = 0.05; 

Periods = 1000./[2 4 8 16 32];

VSdata_old = nan(length(iRS),5,3);
VSdata_new = nan(length(iRS),nStim,3);
Nsp_per_tr = nan(length(iRS),nStim);
FF_stim    = nan(length(iRS),nStim);
for ii = 1:length(iRS)
    
    if strcmp(whichStim,'AC')
        VSdata_old(ii,:,:) = UnitData(theseCells(iRS(ii))).VSdata_spk(:,2:6)';
    end
    
    for ist=1:nStim
        
        % This data
        CSdata = permute(CTTS(iRS(ii),:,:,ist),[3 2 1]);
        
        clear x y ntr VS RS p_VS
        [x,y] = find(CSdata==1);
        
        % Calculate VS from class data
        if strcmp(whichStim,'AC')
            if ist>1 && ist <7
                [VS,RS,p_VS] = vectorstrength(y,Periods(ist-1));
                VSdata_new(ii,ist,:) = [VS RS p_VS];
            end
        end
        
        % Trial similarity measures
        ntr = sum(~isnan(CSdata(:,500)));
        TrSpks = sum(CSdata,2);
        FF_stim(ii,ist) = var(TrSpks,'omitnan')/mean(TrSpks,'omitnan');
        
        % N spikes        
        Nsp_per_tr(ii,ist) = numel(y)/ntr;
    end 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot matrix of Nspks

figure;
set(gcf,'Position',widescreen)

subplot(7,4,(4*3+1):(4*6))
imagesc( Nsp_per_tr(sortCells,sortStim)' )
cmocean('amp')
% caxis([-0.2 dpmaxval])
set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(sortStim))
hcb=colorbar;
hcb.Label.String = 'Nspks';
xlabel('RS cells')


% Normalized
figure;
set(gcf,'Position',widescreen)

subplot(7,4,(4*3+1):(4*6))
imagesc( (Nsp_per_tr(sortCells,sortStim)./repmat(sum(Nsp_per_tr(sortCells,sortStim),2),1,nStim))' )
cmocean('amp')
% caxis([-0.2 dpmaxval])
set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(sortStim))
hcb=colorbar;
hcb.Label.String = 'Nspks';
xlabel('RS cells')


if strcmp(whichStim,'AC')
    
    %% Plot matrix of VS
    
    figure;
    set(gcf,'Position',widescreen)
    
    subplot(7,4,(4*3+1):(4*6))
    imagesc( VSdata_new(sortCells,sortStim,1)' )
    cmocean('algae')
    % caxis([-0.2 dpmaxval])
    set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(sortStim))
    hcb=colorbar;
    hcb.Label.String = 'VS';
    xlabel('RS cells')
    
    
    %% Plot matrix of Rayleigh Stat
    
    figure;
    set(gcf,'Position',widescreen)
    
    subplot(7,4,(4*3+1):(4*6))
    imagesc( VSdata_new(sortCells,sortStim,2)' )
    cmocean('tempo')
    caxis([0 200])
    set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(sortStim))
    hcb=colorbar;
    hcb.Label.String = 'RS';
    xlabel('RS cells')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT -- corr with VS

if strcmp(whichStim,'AC')
    
    VSdata = VSdata_new;
    
    figure;
    set(gcf,'Position',fullscreen)
    for ist=2:6
        
        isp=ist;
        if ist>4
            isp=ist+1;
        end
        subplot(2,5,isp)
        hold on
        
        isig = VSdata(:,ist,3)<alfa;
        
        plot(VSdata(~isig,ist,1),dpStim(~isig,ist),'x','Color',0.6*[1 1 1])
        plot(VSdata(isig,ist,1),dpStim(isig,ist),'ob')
        
        xlabel('VS')
        ylabel(['d'' for ' istlab{ist}])
        ylim([-0.5 4])
        set(gca,'Color','none')
        axis square
        box off
        
        [rs,ps] = corr(VSdata(isig,ist,1),dpStim(isig,ist));
        [ra,pa] = corr(VSdata(:,ist,1),dpStim(:,ist));
        title(sprintf('%s\nsig: r=%0.2f, p=%0.2f\nall: r=%0.2f, p=%0.2f',istlab{ist},rs,ps,ra,pa))
    end
    
    % Mean across stim
    subplot(2,5,10)
    plot(mean(VSdata(:,:,1),2,'omitnan'),mean(dpStim,2),'ok')
    xlabel('mean VS')
    ylabel('mean d''')
    ylim([-0.5 4])
    set(gca,'Color','none')
    axis square
    box off
    [rp,pp] = corr(mean(VSdata(:,:,1),2,'omitnan'),mean(dpStim(:,2:6),2));
    [ra,pa] = corr(mean(VSdata(:,:,1),2,'omitnan'),mean(dpStim,2));
    title(sprintf('Pdc stim: r=%0.2f, p=%0.2e\nAll stim: r=%0.2f, p=%0.2e',rp,pp,ra,pa))
    
    suptitle('Vector strength')
    print_eps_kp(gcf,fullfile(savedir,'SUdp_stim_VScorr'))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT -- corr with RS

if strcmp(whichStim,'AC')
    
    figure;
    set(gcf,'Position',fullscreen)
    
    for ist=2:6
        
        isp=ist;
        if ist>4
            isp=ist+1;
        end
        subplot(2,5,isp)
        hold on
        
        isig = VSdata(:,ist,3)<alfa;
        
        plot(VSdata(~isig,ist,2),dpStim(~isig,ist),'x','Color',0.6*[1 1 1])
        plot(VSdata(isig,ist,2),dpStim(isig,ist),'ob')
        
        xlabel('Rayleigh Statistic')
        ylabel(['d'' for ' istlab{ist}])
        ylim([-0.5 4])
        set(gca,'Color','none')
        axis square
        box off
        
        [rs,ps] = corr(VSdata(isig,ist,2),dpStim(isig,ist));
        [ra,pa] = corr(VSdata(:,ist,2),dpStim(:,ist));
        title(sprintf('%s\nsig: r=%0.2f, p=%0.2f\nall: r=%0.2f, p=%0.2f',istlab{ist},rs,ps,ra,pa))
    end
    
    % Mean across stim
    subplot(2,5,10)
    plot(mean(VSdata(:,:,2),2,'omitnan'),mean(dpStim,2),'ok')
    xlabel('mean RS')
    ylabel('mean d''')
    ylim([-0.5 4])
    set(gca,'Color','none')
    axis square
    box off
    [rp,pp] = corr(mean(VSdata(:,:,2),2,'omitnan'),mean(dpStim(:,2:6),2));
    [ra,pa] = corr(mean(VSdata(:,:,2),2,'omitnan'),mean(dpStim,2));
    title(sprintf('Pdc stim: r=%0.2f, p=%0.2e\nAll stim: r=%0.2f, p=%0.2e',rp,pp,ra,pa))
    
    suptitle('Rayleigh statistic')
    print_eps_kp(gcf,fullfile(savedir,'SUdp_stim_RScorr'))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT -- corr with nspk

figure;
set(gcf,'Position',fullscreen)
for ist=1:nStim
    
    isp=ist;
    if ist>4
        isp=ist+1;
    end
    subplot(2,5,isp)
    hold on
        
    plot(Nsp_per_tr(:,ist),dpStim(:,ist),'ob')
    
    xlabel('Nspk')
    ylabel(['d'' for ' istlab{ist}])
    ylim([-0.5 4])
    set(gca,'Color','none')
    axis square
    box off
    
    [ra,pa] = corr(Nsp_per_tr(:,ist),dpStim(:,ist));
    title(sprintf('%s\nr=%0.2f, p=%0.2e',istlab{ist},ra,pa))
end

% Mean across stim
subplot(2,5,10)
plot(mean(Nsp_per_tr,2),mean(dpStim,2),'ok')
xlabel('mean Nspk')
ylabel('mean d''')
ylim([-0.5 4])
set(gca,'Color','none')
axis square
box off
[rp,pp] = corr(mean(Nsp_per_tr,2),mean(dpStim(:,2:6),2));
[ra,pa] = corr(mean(Nsp_per_tr,2),mean(dpStim,2));
title(sprintf('Pdc stim: r=%0.2f, p=%0.2e\nAll stim: r=%0.2f, p=%0.2e',rp,pp,ra,pa))

suptitle('N spikes')
print_eps_kp(gcf,fullfile(savedir,'SUdp_stim_FRcorr'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cumulative fraction of spikes -- cell/stim pairs

[dpALLsort,iALLsort] = sort(dpStim(:),'ascend');
NspkALL = Nsp_per_tr(:);

figure; hold on
yyaxis left
plot(1:numel(NspkALL),cumsum(NspkALL(iALLsort))./sum(NspkALL(iALLsort)),'b','LineWidth',3)
box off
set(gca,'Color','none')
yyaxis right
plot(1:numel(NspkALL),dpALLsort,'.k','MarkerSize',15)
xlim([1 numel(NspkALL)])
box off
set(gca,'Color','none')

title('Cumulative fraction of spikes -- cell/stim pairs')
print_eps_kp(gcf,fullfile(savedir,'SU_cumulNspks_CSpairs'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT -- corr with FF

figure;
set(gcf,'Position',fullscreen)
for ist=1:nStim
    
    isp=ist;
    if ist>4
        isp=ist+1;
    end
    subplot(2,5,isp)
    hold on
        
    plot(FF_stim(:,ist),dpStim(:,ist),'ob')
    
    xlabel('FF')
    ylabel(['d'' for ' istlab{ist}])
    ylim([-0.5 4])
    set(gca,'Color','none')
    axis square
    box off
    
    [ra,pa] = corr(FF_stim(~isnan(FF_stim(:,ist)),ist),dpStim(~isnan(FF_stim(:,ist)),ist));
    title(sprintf('%s\nr=%0.2f, p=%0.2f',istlab{ist},ra,pa))
end

% Mean across stim
subplot(2,5,10)
plot(mean(FF_stim,2,'omitnan'),mean(dpStim,2),'ok')
xlabel('mean FF')
ylabel('mean d''')
ylim([-0.5 4])
set(gca,'Color','none')
axis square
box off
[rp,pp] = corr(mean(FF_stim,2,'omitnan'),mean(dpStim(:,2:6),2));
[ra,pa] = corr(mean(FF_stim,2,'omitnan'),mean(dpStim,2));
title(sprintf('Pdc stim: r=%0.2f, p=%0.2f\nAll stim: r=%0.2f, p=%0.2f',rp,pp,ra,pa))

suptitle('Fano factor')
print_eps_kp(gcf,fullfile(savedir,'SUdp_stim_FFcorr'))












