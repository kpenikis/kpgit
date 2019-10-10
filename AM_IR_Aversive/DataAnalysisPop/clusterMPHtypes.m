function clusterMPHtypes
% 
% clusterMPHtypes
%
%  Plot MPHs, and cluster by response shape.
% 
%  Intended to help categorize response types.
% 

global AMrates 

close all
    
%% Load saved MPH data

fn = set_paths_directories('','',1);
loaddir = fullfile(fn.figs,'PopulationTuning');
if ~exist(loaddir,'dir')
    mkdir(loaddir)
end

% load(fullfile(loaddir,'MPHdata_ownSpkShifts')); %_ownSpkShifts
load(fullfile(loaddir,'MPHdata')); 


%% Load Unit data files

q = load(fullfile(fn.processed,'Units_250'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load RCorr results
% q=load(fullfile(fn.figs,'RCorr','exclOnset','PCMat_10trTemp_b.mat'));
% PCMat = q.PCMat;
% clear q


%% Set up

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',8)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

keyboard


%% Compare latencies with common or unique spike shifts

qOS=load(fullfile(loaddir,'MPHdata_ownSpkShifts')); %_ownSpkShifts
qCS=load(fullfile(loaddir,'MPHdata')); 

clipVal = 2;

for ir = 1:5
    
    iUn_GoodData = find(sum(~isnan(qOS.zFR_vec(:,:,ir)),2)==ceil(1000/AMrates(ir)));
    
% % %     % Filter to just lpn=5000, spl=45
% % %     iUn_GoodData = iUn_GoodData([UnitData(iUn_GoodData).lpn]==5000 & [UnitData(iUn_GoodData).spl]==60);
    
    % Own spike shift values
    % Remove empty units, clip at limit
    GoodData = qOS.zFR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);
    GoodData(GoodData> clipVal)  =  clipVal;
    GoodData(GoodData<-clipVal)  = -clipVal;
    
    % Same manipulations for common shift data
    ComSS_Data = qCS.zFR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);
    ComSS_Data(ComSS_Data> clipVal)  =  clipVal;
    ComSS_Data(ComSS_Data<-clipVal)  = -clipVal;
    
    % Label
    flagNS = UnitInfo.TroughPeak(iUn_GoodData)<0.5;
    
    % Calulate linkages
    Y = pdist(GoodData);     % cells X time
    Z = linkage(Y,'ward');
    
    leafOrder = optimalleaforder(Z,Y);
    
    % Plot dendrogram
    figure;
    set(gcf,'Position',fullscreen)
    
    subplot(1,6,3:6);      % ,'ColorThreshold',40 ,'reorder',leafOrder
    [hd,tvals,outperm] = dendrogram(Z,0,'reorder',leafOrder,'ColorThreshold',10,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
    set(gca,'tickdir','out')
    
    RespCluLabels = cluster(Z,'maxclust',10);
    
    % Plot MPH - OWN spike shift
    subplot(1,6,2);
    
    imagesc( GoodData(fliplr(outperm),:) )
    
    caxis([-clipVal clipVal])
    cmocean('balance','pivot',0)
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0 size(GoodData,1)+1])
    set(gca,'ytick',find(flagNS(fliplr(outperm))),'yticklabel','N','xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    title('Own spike shift')
    
    % Plot MPH - COMMON spike shift
    subplot(1,6,1);
    
    imagesc( ComSS_Data(fliplr(outperm),:) )
    
    caxis([-clipVal clipVal])
    cmocean('balance','pivot',0)
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0 size(GoodData,1)+1])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    title('Common spike shift')
    
    suptitle([num2str(AMrates(ir)) ' Hz'])
    
% %     % Save
    savedir = fullfile(loaddir,'Pop_MPH');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(gcf,fullfile(savedir,['CompareShifts_MPH_' num2str(AMrates(ir)) 'Hz_clip' num2str(clipVal)]))
    
end


%% Cluster response shapes  ---  pca 

for ir = 1:5
    
    mph_pca
    
end


%% Cluster response shapes  ---  raw zFR

for ir=1:5
%     T = clusterdata(zFR_vec(:,:,ir),10);
        
    % Remove empty units
    iUn_GoodData = find(sum(~isnan(zFR_vec(:,:,ir)),2)==ceil(1000/AMrates(ir)));
    GoodData = zFR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);
    
    % Calulate linkages
    Y = pdist(GoodData);     % cells X time
    Z = linkage(Y,'ward');
    
    leafOrder = optimalleaforder(Z,Y);
    
    % Plot dendrogram
    figure;
    set(gcf,'Position',fullscreen)
    
    subplot(1,6,2:6)
    [hd,tvals,outperm] = dendrogram(Z,0,'reorder',leafOrder,'ColorThreshold',10,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
    set(gca,'tickdir','out')
    
    RespCluLabels = cluster(Z,'maxclust',10);
    
    % Plot MPH
    subplot(1,6,1);
    
    imagesc( GoodData(fliplr(outperm),:) )
    
    caxis([-2 6])
    cmocean('balance','pivot',0)
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0 size(GoodData,1)+1])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    
    suptitle([num2str(AMrates(ir)) ' Hz'])
end


%% Cluster response shapes  ---  clipped zFR

clipVal = 1;

for ir = 1:5
%     T = clusterdata(zFR_vec(:,:,ir),10);
        
    % Remove empty units
    iUn_GoodData = find(sum(~isnan(zFR_vec(:,:,ir)),2)==ceil(1000/AMrates(ir)));
    GoodData = zFR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);
    GoodData(GoodData> clipVal)  =  clipVal;
    GoodData(GoodData<-clipVal)  = -clipVal;
    
    %%% DOUBLE THE PERIOD
    GoodData = [GoodData GoodData];
    
    % Calulate linkages
    Y = pdist(GoodData);     % cells X time
    Z = linkage(Y,'ward');
    
    leafOrder = optimalleaforder(Z,Y);
    
    % Plot dendrogram
    figure;
    set(gcf,'Position',fullscreen)
    
    subplot(1,6,2:6);      % ,'ColorThreshold',40 ,'reorder',leafOrder
    [hd,tvals,outperm] = dendrogram(Z,0,'reorder',leafOrder,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
    set(gca,'tickdir','out')
    
    RespCluLabels = cluster(Z,'maxclust',10);
    
    % Plot MPH
    subplot(1,6,1);
    
    imagesc( GoodData(fliplr(outperm),:) )
    
    caxis([-clipVal clipVal])
    cmocean('balance','pivot',0)
    xlim([0 ceil(1000/AMrates(ir))*2])
    ylim([0 size(GoodData,1)+1])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    
    suptitle([num2str(AMrates(ir)) ' Hz'])
    
    % Save
    savedir = fullfile(loaddir,'Pop_MPH');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(gcf,fullfile(savedir,['MPH_double_' num2str(AMrates(ir)) 'Hz_clip1_commonShift']))
    
end


%% Cluster response shapes  ---  zFR normalized each un

for ir=1:5
%     T = clusterdata(zFR_vec(:,:,ir),10);
    
    % Remove empty units
    iUn_GoodData = find(sum(~isnan(FR_vec(:,:,ir)),2)==ceil(1000/AMrates(ir)));
    GoodData = FR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);
    
    % Option: normalize to range [0 1]
    foo =  (GoodData-min(GoodData,[],2)) ;
    NormData =  foo ./ (max(GoodData,[],2) + 0.1);
    
    % Option: log transform
    GoodData = log10(GoodData+0.1);
    
    % Calulate linkages
    Y = pdist(GoodData);     % cells X time
    Z = linkage(Y,'ward');
    
    leafOrder = optimalleaforder(Z,Y);
    
    % Plot dendrogram
    figure;
    set(gcf,'Position',fullscreen)
    
    subplot(1,6,2:6)    % ,'reorder',leafOrder
    [hd,tvals,outperm] = dendrogram(Z,0,'ColorThreshold',50,'reorder',leafOrder,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
    set(gca,'tickdir','out')
    
    RespCluLabels = cluster(Z,'maxclust',10);
    
    % Plot MPH
    subplot(1,6,1);
    
    imagesc( GoodData(fliplr(outperm),:) )
    
    cmocean('gray')
    caxis([-1 1.75])
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0 size(GoodData,1)+1])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    
    suptitle([num2str(AMrates(ir)) ' Hz'])
    
% % %     % Save
%     savedir = fullfile(loaddir,'Pop_MPH');
%     if ~exist(savedir,'dir')
%         mkdir(savedir)
%     end
%     print_eps_kp(gcf,fullfile(savedir,['logFR_MPH_' num2str(AMrates(ir)) 'Hz']))
    
end




end





