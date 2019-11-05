function ClassStim_plotResults
% ClassStim_plotResults
% 
%  Load and plot stimulus classification results.
% 
%  KP, 2019-10
%


close all
global fn trMin nIterations StimDur exclOnset TempSize

%~~~~~~~~~~~~~~~~~~
fl_SUDist   = 0;
%~~~~~~~~~~~~~~~~~~
fl_PCmat    = 0;
%~~~~~~~~~~~~~~~~~~
fl_FF_Acc   = 0;
%~~~~~~~~~~~~~~~~~~
fl_PCdiffs  = 0;
%~~~~~~~~~~~~~~~~~~
fl_stimHist = 1;
%~~~~~~~~~~~~~~~~~~

exclOnset   = 0; 
StimDur     = 500;
if exclOnset>0
    StimDur = StimDur-exclOnset;
end
app_str   = '_15trTemp';


%% Load data

fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'StimClass');
if exclOnset>0
    savedir = fullfile(savedir,'exclOnset');
end

q=load(fullfile(savedir,['PCMat' app_str]));
PCMat  = q.PCMat;

% q=load(fullfile(savedir,['PCMat_tsPop_rcorr' app_str '_' num2str(StimDur) 'ms']));
% PCMat_Pop = q.PCMat;

% Get one Info struct
load(fullfile(fn.processed,'AAB_265054',sprintf('%s_sess-%s_Info','AAB_265054','Apr02-AM')))


%% Plot results

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
quartscreen = [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];


%% ~~~~ Distribution of PC vals for each stimulus 

if fl_SUDist 
    
    hf1 = figure;
    set(hf1,'Position',quartscreen,'NextPlot','add')
    plot([0 size(PCMat,1)+1],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
%     plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.8 0],'LineWidth',2)
    hold on
%     plot([0 10],[1 1]./size(PCMat_Pop,1),'Color',[0 1 0.8],'LineWidth',2)
    
    PCs=[];
    for stid = 1:size(PCMat,1)
        
        for iUn = 1:size(PCMat,3)
            PCs = [PCs; stid PCMat(stid,stid,iUn)];
        end
        
        % Manually make boxplot
        q5  = quantile(PCs(PCs(:,1)==stid,2),0.05);
        q25 = quantile(PCs(PCs(:,1)==stid,2),0.25);
        q75 = quantile(PCs(PCs(:,1)==stid,2),0.75);
        q95 = quantile(PCs(PCs(:,1)==stid,2),0.95);
        
        plot(stid*[1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
        fill(stid+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')
        
    end
    
    plotSpread(PCs(:,2),'distributionIdx',PCs(:,1),'distributionColors','k','showMM',3)
    
    % Plot population nSpk vector results
%     plot(diag(PCMat_Pop)','o','Color',[0 1 0.8],'MarkerSize',12,'LineWidth',4)
    % plot(diag(PCMat_Vec)','o','Color',[0 1 0.8],'MarkerSize',12,'LineWidth',4)
    
    xlim([0 size(PCMat,1)+1])
    ylim([0 1])
    ylabel('Percent Assigned Correctly')
    set(gca,'xtick',1:size(PCMat,1),'xticklabel',Info.stim_ID_key',...
        'Color','none','tickdir','out')
    title(['Stimulus classification, SU data' app_str])
    
    % Save figure
    print_eps_kp(hf1,fullfile(savedir,['PC_SUdistr' app_str]))
    
end


%% ~~~~ Average PC mat

if fl_PCmat
    
    hf = figure;
    imagesc(mean(PCMat,3,'omitnan'))
    axis square
    caxis([0 1])
    cmocean('ice')
    colorbar
    ylabel('True stim')
    xlabel('Assigned')
    title(['PC_meanMat' app_str ''])
    
    print(fullfile(savedir,['PC_meanMat' app_str]),'-dpdf')
    
end


%% ~~~~ Relationship of PC to FF
if fl_FF_Acc 
    
    load(fullfile(fn.figs,'FF','FF_Stim1-8_0-500ms'))
    
    hfFF=figure;
    hold on 
    set(gcf,'Position',fullscreen)
    
    % Each datapoint: average per unit across stim
    meanAccuracyUn = nan(1,size(PCMat,3));
    for iUn = 1:size(PCMat,3)
        meanAccuracyUn(iUn) = mean(diag(PCMat(:,:,iUn)),'omitnan');
    end
    meanFanoFactUn = mean(FanoFs,2,'omitnan')';
    
    subplot(2,3,1); 
    plot([0 5],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
    hold on
    plot(meanFanoFactUn,meanAccuracyUn,'ok')
    xlabel('Fano factor')
    ylabel('Stim Class Accuracy')
    ylim([0 1])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title('Average for each unit')
    % Stats
    safe_dps = sum([isnan(meanFanoFactUn)' isnan(meanAccuracyUn)'],2)==0;
    [rP,pP] = corr(meanAccuracyUn(safe_dps)',meanFanoFactUn(safe_dps)','type','Pearson');
    [rS,pS] = corr(meanAccuracyUn(safe_dps)',meanFanoFactUn(safe_dps)','type','Spearman');
    text(3,0.5,sprintf('Pearson: r= %0.2f, p=%0.2e\nSpearman: r= %0.2f, p=%0.2e',rP,pP,rS,pS))
    
    
    % Each datapoint: one unit/stim
    FFData  = [];
    AccData = [];
    for ist = 1:size(PCMat,1)
        FFData  = [FFData; FanoFs(:,ist)];
        AccData = [AccData; permute(PCMat(ist,ist,:),[3,1,2])];
    end
    
    subplot(2,3,2); 
    plot([0 10],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
    hold on
    plot(FFData,AccData,'ok')
    xlabel('Fano factor')
    ylabel('Stim Class Accuracy')
    ylim([0 1])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title('Each stim for each unit')
    % Stats
    safe_dps = sum([isnan(FFData) isnan(AccData)],2)==0;
    [rP,pP] = corr(FFData(safe_dps),AccData(safe_dps),'type','Pearson');
    [rS,pS] = corr(FFData(safe_dps),AccData(safe_dps),'type','Spearman');
    text(6,0.65,sprintf('Pearson: r= %0.2f, p=%0.2e\nSpearman: r= %0.2f, p=%0.2e',rP,pP,rS,pS))
    
    
    % Each datapoint: one unit/stim
    FFData  = [];
    AccData = [];
    for iUn = 1:size(PCMat,3)
        [mPC,mst] = max(diag(PCMat(:,:,iUn)));
        FFData  = [FFData; FanoFs(iUn,mst)];
        AccData = [AccData; mPC];
    end
    
    subplot(2,3,3); 
    plot([0 5],[1 1]./size(PCMat,1),'Color',[1 0.5 0.7],'LineWidth',2)
    hold on
    plot(FFData,AccData,'ok')
    xlabel('Fano factor')
    ylabel('Stim Class Accuracy')
    ylim([0 1])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title('Best stim for each unit')
    % Stats
    safe_dps = sum([isnan(FFData) isnan(AccData)],2)==0;
    [rP,pP] = corr(FFData(safe_dps),AccData(safe_dps),'type','Pearson');
    [rS,pS] = corr(FFData(safe_dps),AccData(safe_dps),'type','Spearman');
    text(3,0.7,sprintf('Pearson: r= %0.2f, p=%0.2e\nSpearman: r= %0.2f, p=%0.2e',rP,pP,rS,pS))
    
    % Save figure
    print_eps_kp(hfFF,fullfile(savedir,['PC_vsFF' app_str]))
    
end




%% DIFFERENCES

if fl_PCdiffs
    
    load(fullfile(fn.figs,'FF','FF_Stim1-8_0-500ms'))
    
q=load(fullfile(savedir,['PCMat' '_15trTemp']));
PCMat_15  = q.PCMat;
q=load(fullfile(savedir,['PCMat' '_1trTemp']));
PCMat_1   = q.PCMat;
clear q

PCdiffs     = nan(size(PCMat_15,3),size(PCMat_15,1));
PCnormdiffs = nan(size(PCMat_15,3),size(PCMat_15,1));
for iUn = 1:size(PCMat_15,3)
    
    PCdiffs(iUn,:)     = [diag(PCMat_15(:,:,iUn)) - diag(PCMat_1(:,:,iUn))]';
    PCnormdiffs(iUn,:) = [(diag(PCMat_15(:,:,iUn)) - diag(PCMat_1(:,:,iUn)))./diag(PCMat_1(:,:,iUn))]';
    
end

% All are about the same; best combo of the following (by eye) is median FF vs mean raw diff

hfda = figure;
plot(FanoFs(:),PCnormdiffs(:),'ok')
xlabel('Fano factor')
ylabel('PC 15tr - PC 1tr')
% ylim([0 1])
set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
title('All stimuli for all units')


xdata = median(FanoFs,2,'omitnan');
ydata = median(PCnormdiffs,2,'omitnan');

hfdm = figure;
plot(xdata,ydata,'ok')
xlabel('Fano factor')
ylabel('PC 15tr - PC 1tr')
% ylim([0 1])
set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
title('Average for each unit')
% Stats
safe_dps = sum([isnan(xdata) isnan(ydata)],2)==0;
[rP,pP] = corr(xdata(safe_dps),ydata(safe_dps),'type','Pearson');
[rS,pS] = corr(xdata(safe_dps),ydata(safe_dps),'type','Spearman');
text(2,0.28,sprintf('Pearson: r= %0.2f, p=%0.2e\nSpearman: r= %0.2f, p=%0.2e',rP,pP,rS,pS))
    

% PA matrix
hmd=figure;
imagesc(mean(PCMat_15,3,'omitnan')-mean(PCMat_1,3,'omitnan'))
% imagesc((mean(PCMat_15,3,'omitnan')-mean(PCMat_1,3,'omitnan'))./mean(PCMat_1,3,'omitnan').*100)
axis square
caxis(0.2*[-1 1])
cmocean('curl','pivot',0) 
colorbar('TickDirection','out')
ylabel('True stim')
xlabel('Assigned')
title('PA difference, 15tr-1tr')

% Save figure
print_eps_kp(hmd,fullfile(savedir,'PC_diffMat'))
print(fullfile(savedir,'PC_diffMat'),'-dpdf')


end




%% # stim > thresh histogram

if fl_stimHist
    
q=load(fullfile(savedir,['PCMat' '_15trTemp']));
PCMat_15  = q.PCMat;
q=load(fullfile(savedir,['PCMat' '_1trTemp']));
PCMat_1   = q.PCMat;
clear q


Thresholds = 4/size(PCMat_1,1);%0.2:0.1:0.5;

for iTh = 1:numel(Thresholds)
    
    AccuracyThreshold = Thresholds(iTh);
    
    nStim_Classified_15 = nan(size(PCMat_15,3),1);
    nStim_Classified_1  = nan(size(PCMat_1 ,3),1);
    
    for iUn = 1:size(PCMat_1,3)
        nStim_Classified_15(iUn) = sum(diag(PCMat_15(:,:,iUn))>AccuracyThreshold);
        nStim_Classified_1(iUn)  = sum(diag(PCMat_1(:,:,iUn))>AccuracyThreshold);
    end
    
    hf(iTh)=figure;
    subplot(2,1,1);
    histogram(nStim_Classified_1,'FaceColor',[1 1 1].*0.7)
    hold on
    text(-0.2,65,num2str(sum(nStim_Classified_1==0)))
    ylim([0 70])
    xlim([-1 9])
    set(gca,'xtick',0:8,'tickdir','out','ytick',0:10:70,'Color','none')
    title('1 tr template')
    
    subplot(2,1,2);
    histogram(nStim_Classified_15,'FaceColor',[1 1 1].*0.7)
    hold on
    text(-0.2,65,num2str(sum(nStim_Classified_15==0)))
    ylim([0 70])
    xlim([-1 9])
    set(gca,'xtick',0:8,'tickdir','out','ytick',0:10:70,'Color','none')
    title('15 tr template')
    
    suptitle(['Number of stim classfied above ' num2str(AccuracyThreshold) ' PC'])
    
    print_eps_kp(hf(iTh), fullfile(savedir,sprintf('nStimHist_thresh%i',100*AccuracyThreshold)))
    
end


end %fl_stimHist






end