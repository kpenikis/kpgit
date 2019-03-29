function inspectTrialActivity(SUBJECT,SESSION)
% 
% Must be single session.
% Gathers raster data for all units and all stimuli. 
% Then look at them in a bunch of different ways.
% 
% KP, 2019-03
% 


scrsz = get(0,'ScreenSize');   %[left bottom width height]
widescreen  = [1 scrsz(4) scrsz(3) scrsz(4)/2];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)/2];

close all
fn = set_paths_directories(SUBJECT,SESSION,1);

% Get raster data
[UnitInfo, UnitData, Info, TrialData, Clusters, StimResp ] = collectRasterDataSession(SUBJECT,SESSION);



%% Look at trial-by-trial data in various ways

keyboard 

% Select stimulus to inspect
stid = 3;

% Animated view of responses per trial across all units
% hf=figure;
% for it = 1:size(StimResp(stid).raster,1)
%     imagesc( flipud(permute(StimResp(stid).raster(it,:,:),[3,2,1])) )
%     %         c=flipud(gray);
%     colormap(gray)
%     pause(2)
% end
% close(hf);


% Looking for relationship between total number of spikes per trial across
% units
nSpksTr = permute( sum(StimResp(stid).raster,2), [3,1,2]);
nSpksTr = flipud(nSpksTr);

% Plot nspikes on each trial
figure; 
plot(nSpksTr,'LineWidth',2)
set(gca,'yscale','log')
xlabel('Units')
ylabel('nSpks/trial')

% Plot trials ranked by n spikes for each cell
trRank = zeros(size(nSpksTr));
for ii=1:size(nSpksTr,1)
    [~,trRank(ii,:)] = sort(nSpksTr(ii,:),'descend');
end
trRank = flipud(trRank);

figure;
plot(trRank)
xlabel('Units')
ylabel('trial rank')


% Select a unit to sort n spikes per trial, and plot all
[~,idx] = sort(nSpksTr(3,:),'descend');
figure; 
ic=imagesc(log(nSpksTr(:,idx)));
colormap(gray)
xlabel('Trials, ordered by nspks of selected unit')
ylabel('Units, sorted by baseline FR')
title('Total number of spikes per trial for each unit')

% Same, but RANK not nspikes  -- some pattern here
[~,idx] = sort(trRank(3,:),'ascend');
figure; 
ic=imagesc(trRank(:,idx));
colormap(gray)
xlabel('Trials, ordered by nspk ranked trs of selected unit')
ylabel('Units, sorted by baseline FR')
title('Trials ranked by nspks for each unit')


% Pick 2 units to compare  --  no correlation of nspikes across trials for many units
un1 = 3; un2 = 4;
figure;
scatter(nSpksTr(un1,:),nSpksTr(un2,:),'k')
axis square
ymaxval = 50;
set(gca,'xlim',[0 ymaxval],'ylim',[0 ymaxval])
[r,p]=corrcoef(nSpksTr(un1,:),nSpksTr(un2,:));


% Compare sum of NS and RS n spikes (correlated for 4Hz, not correlated for 2+8Hz + AC)
figure;
scatter(sum(nSpksTr(1:2,:),1),sum(nSpksTr(3:end,:)),'k')
axis square
ymaxval = 110;
set(gca,'xlim',[0 ymaxval],'ylim',[0 ymaxval])
xlabel('sum NS spikes per trial')
ylabel('sum RS spikes per trial')
[r,p]=corrcoef(sum(nSpksTr(1:2,:),1),sum(nSpksTr(3:end,:)));


%  
[r,p]=corr(nSpksTr');
ins = find(p>0.05);
r(ins)=0;

figure;
imagesc(r)
axis square
colormap(cmocean('curl'))
colorbar



%% Average tr/tr correlation of responses to each stimulus
%  within ONE unit

iUn = size(UnitInfo,1);

for stid = unique([TrialData.trID])'
    
    if stid==0, continue, end
    
    thisRaster = UnitData(iUn).raster{stid}';
    
    [r,p] = corrcoef(thisRaster);
    
    idx = eye(size(r));
    idx = ~idx;
    
    r_idx=r(idx);
    ProportionSigCorrTrs = sum(p(idx)<0.05)/numel(p(idx));
    DistrSigTrCorrs = r_idx(p(idx)<0.05);
    
    % Shuffled spike trains, for control
    shuffledRaster = nan(size(thisRaster));
    for it = 1:size(thisRaster,2)
        shuffledRaster(:,it) = thisRaster(randperm(size(thisRaster,1)),it);
    end
    
    [r,p] = corrcoef(shuffledRaster);
    
    r_idx=r(idx);
    ProportionSigCorrTrs_shuff = sum(p(idx)<0.05)/numel(p(idx));
    DistrSigTrCorrs_shuff = r_idx(p(idx)<0.05);
    
    hf(stid)=figure;
    set(hf(stid),'Position',widescreen,'NextPlot','add')
    
    subplot(1,4,1); hold on
    ib(1)=bar(1,ProportionSigCorrTrs,'FaceColor','b','EdgeColor','none');
    ib(2)=bar(2,ProportionSigCorrTrs_shuff,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
    set(gca,'xtick',[1 2],'xticklabels',{'Data' 'Shuffled'},...
        'xlim',[0.5 2.5],'ylim',[0 0.2])
    ylabel(sprintf('Proportion of trial comparisons\nwith spiketrain correlation p<0.05'))
    
    subplot(1,4,2:4); hold on
    ih(1) = histogram(DistrSigTrCorrs,'BinEdges',0:0.02:0.3,'FaceColor','b','EdgeColor','none');
    ih(2) = histogram(DistrSigTrCorrs_shuff,'BinEdges',0:0.02:0.3,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
    xlim([0 0.3])
    ylabel('Count')
    xlabel('Significant correlation values')
    
    suptitle(sprintf('%s stimulus, Clu #%i',Info.stim_ID_key{stid},UnitData(iUn).Clu))
    
    
end %stid


%% Average correlation of each trial to the mean response of each stimulus
% first for each unit

% Preallocate: each Un
ProportionSigCorrTrs = nan( numel(UnitData), numel(unique([TrialData.trID]))-1 );
MeanCorrSigTr        = nan( numel(UnitData), numel(unique([TrialData.trID]))-1 );
STDCorrSigTr         = nan( numel(UnitData), numel(unique([TrialData.trID]))-1 );

% Preallocate: aggregate activity
PSigCorrTrs_all      = nan( 1, numel(unique([TrialData.trID]))-1 );
MeanCorrSigTr_all    = nan( 1, numel(unique([TrialData.trID]))-1 );
STDCorrSigTr_all     = nan( 1, numel(unique([TrialData.trID]))-1 );

% Preallocate: NS aggregate
PSigCorrTrs_NS       = nan( 1, numel(unique([TrialData.trID]))-1 );
MeanCorrSigTr_NS     = nan( 1, numel(unique([TrialData.trID]))-1 );
STDCorrSigTr_NS      = nan( 1, numel(unique([TrialData.trID]))-1 );

% Preallocate: RS aggregate
PSigCorrTrs_RS       = nan( 1, numel(unique([TrialData.trID]))-1 );
MeanCorrSigTr_RS     = nan( 1, numel(unique([TrialData.trID]))-1 );
STDCorrSigTr_RS      = nan( 1, numel(unique([TrialData.trID]))-1 );

hf(stid)=figure;
plot([0 1],[0 1],'-k')
hold on
for stid = unique([TrialData.trID])'
    
    if stid==0, continue, end
    
    for iUn = 1:numel(UnitData)
        
        thisRaster = StimResp(stid).raster(:,:,iUn);
        meanResp   = mean(thisRaster,1);
        
        [r,p] = corr(thisRaster',meanResp');
        
        ProportionSigCorrTrs(iUn,stid) = sum(p<0.05)/numel(p);
        MeanCorrSigTr(iUn,stid)        = mean(r(p<0.05));
        STDCorrSigTr(iUn,stid)         = std(r(p<0.05));
        
    end %iUn
    
    % Aggregate activity of all units
    thisRaster = sum(StimResp(stid).raster,3);
    meanResp   = mean(thisRaster,1);
    
    [r,p] = corr(thisRaster',meanResp');
    
    PSigCorrTrs_all(stid)     = sum(p<0.05)/numel(p);
    MeanCorrSigTr_all(stid)   = mean(r(p<0.05));
    STDCorrSigTr_all(stid)    = std(r(p<0.05));
    
    % Aggregate activity of NS units
    thisRaster = sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3);
    meanResp   = mean(thisRaster,1);
    
    [r_NS,p] = corr(thisRaster',meanResp');
    
    PSigCorrTrs_NS(stid)     = sum(p<0.05)/numel(p);
    MeanCorrSigTr_NS(stid)   = mean(r_NS(p<0.05));
    STDCorrSigTr_NS(stid)    = std(r_NS(p<0.05));
    
    % Aggregate activity of RS units
    thisRaster = sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3);
    meanResp   = mean(thisRaster,1);
    
    [r_RS,p] = corr(thisRaster',meanResp');
    
    PSigCorrTrs_RS(stid)     = sum(p<0.05)/numel(p);
    MeanCorrSigTr_RS(stid)   = mean(r_RS(p<0.05));
    STDCorrSigTr_RS(stid)    = std(r_RS(p<0.05));
    
    plot(r_NS,r_RS,'ob','MarkerSize',10)
    axis square
    set(gca,'xlim',[0 0.5],'ylim',[0 0.5])
    xlabel('Avg correlation of NS units this trial')
    ylabel('Avg correlation of RS units')
    
end %stid

hf2=figure;
plot(repmat(1:7,size(MeanCorrSigTr,1),1)', MeanCorrSigTr',...
    'ok','LineWidth',1,'MarkerSize',5)
hold on
plot((1:7) + 0.15,MeanCorrSigTr_NS,'xr','MarkerSize',10,'LineWidth',2)
plot((1:7) + 0.15,MeanCorrSigTr_RS,'xb','MarkerSize',10,'LineWidth',2)
plot((1:7) + 0.3,MeanCorrSigTr_all,'.g','MarkerSize',30)
set(gca,'xlim',[0 8],'ylim',[0 1],'xtick',1:7,'xticklabel',Info.stim_ID_key(1:7))
ylabel('Average correlation of each trial to the mean activity')





%% Fano Factor

FF     = nan( numel(UnitData), numel(unique([TrialData.trID])) );
FF_all = nan( 1, numel(unique([TrialData.trID])) );
FF_NS  = nan( 1, numel(unique([TrialData.trID])) );
FF_RS  = nan( 1, numel(unique([TrialData.trID])) );

for stid = unique([TrialData.trID])'
    if stid==0
        FF(:,1)     =[];
        FF_all(:,1) =[];
        FF_NS(:,1)  =[];
        FF_RS(:,1)  =[];
        continue
    end
    
    for iUn = 1:numel(UnitData)
        FF(iUn,stid) = var(sum(UnitData(iUn).raster{stid},2))/mean(sum(UnitData(iUn).raster{stid},2));
    end
    
    
    % Aggregate population
    FF_all(1,stid) = var( sum( sum(StimResp(stid).raster,3) ,2) ) / mean( sum( sum(StimResp(stid).raster,3) ,2) );
    
    % Just NS
    FF_NS(1,stid) = var(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3),2))/mean(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3),2));
    
    % Just RS
    FF_RS(1,stid) = var(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3),2))/mean(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3),2));
    
    
    % Plot: difference in spikes from mean for stimulus
    hf(stid)=figure;
    plot([-40 80],[-40 80],'k')
    hold on
    plot([0 0],[-40 80],'k')
    plot([-40 80],[0 0],'k')
%     plot(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3),2),...
%         sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3),2),...
%         'ob','MarkerSize',10)
    plot(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3),2) - mean(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'NS')),3),2)),...
        sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3),2) - mean(sum(sum(StimResp(stid).raster(:,:,strcmp(UnitInfo.SpkShape,'RS')),3),2)),...
        'ob','MarkerSize',10)
    axis square
    xlabel('NS')
    ylabel('RS')
    title(Info.stim_ID_key{stid})
    
end %stid


hf1=figure;
bar([nanmean(FF,1)' FF_all' FF_NS' FF_RS'])

hf2=figure;
errorbar(nanmean(FF,1),nanstd(FF,1)./nanmean(FF,1))
set(gca,'xlim',[0 8],'ylim',[0 5])
title('Fano-factor: mean across cells')
set(gca,'xtick',1:7,'xticklabel',Info.stim_ID_key(1:7))



end


