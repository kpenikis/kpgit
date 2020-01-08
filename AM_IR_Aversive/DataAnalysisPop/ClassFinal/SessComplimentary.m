
% close all
whichStim = 'AC'; %finish generalizing

crAmt = 0.0001;

% Data settings
fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'ClassAM',whichStim,'Full','SessSubpops');

Sessions = {...
    'Mid10RS_Apr02';...
    'Mid10RS_Apr11';...
    'Mid10RS_Mar26';...
    'Mid10RS_Mar28';...
    'Mid10RS_Mar30';...
    'Mid10RS_Jan17';...
    'Mid10RS_Jan21';...
    };

% Load SU data
datadir = fullfile(fn.figs,'ClassAM','AC','Full','each');
tabsavename = sprintf('CR_%s.mat','each');
q=load(fullfile(datadir,tabsavename));
CR_SU = q.CR;
clear q

% Load CTTS data
[CTTS,theseCells,nUns,Dur,nStim,TrainSize,TestSize,UnitData] = recallDataParams(whichStim,'each');


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');


%% By session, load data 


for is = 4:numel(Sessions)
    
    % Get Subpop data
    datadir = fullfile(fn.figs,'ClassAM','AC','Full',Sessions{is});
    tabsavename = sprintf('CR_vFull_%s.mat',Sessions{is});
    q=load(fullfile(datadir,tabsavename));
    if size(q.CR,1)>1
        keyboard
    end
    
    ConfMat = mean(q.CR.Results{:},3);
    dpSP    = dp_from_ConfMat(ConfMat,crAmt);
    
    
    % Get SU data
    UnitIdx = q.CR.SUids{:}; %entries of CTTS/theseCells
    
    dpSU = nan(numel(UnitIdx),size(CTTS,4));
    for iu = 1:numel(UnitIdx)
        ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
        dpSU(iu,:) = dp_from_ConfMat(ConfMat,crAmt);
    end
    
    % Get dot product values
    PSTHSize = TrainSize-TestSize;
    [CorrVals,TheseTrials] = checkCorrAct(CTTS(UnitIdx,:,:,:),1,size(CTTS,4),numel(UnitIdx),'sim',PSTHSize);
    
    % Calculate PC for each trial of stim N
    ST = 6;
    TrRes=q.CR.TrRes{:}((q.CR.TrRes{:}(:,1)==ST),:);
    
    StTrs = unique(TrRes(:,3));
    TrPC  = nan(max(StTrs),1);
    for it = StTrs'
        idx = TrRes(:,3)==it;
        TrPC(it) = sum(TrRes(idx,2)==ST)/sum(TrRes(idx,1)==ST);
    end
    
    % Chronological order
    figure;
    
    subplot(4,1,1)
    plot(TrPC,'-k')
    xlim([0.5 size(TrPC,1)+0.5])
    set(gca,'xtick',[])
    title('stim 6')
    
    subplot(4,1,2:4)
    imagesc(CorrVals{ST})
    colormap('gray')
%     colorbar
    set(gca,'clim',[0 1])
    
    
    
    % Sort trials by PC
    TrPC_matched = TrPC(TheseTrials{ST});
    [tPC,itPC] = sort(TrPC_matched,'descend');
    
%     [tPC,itPC] = sort(TrPC,'descend');
%     itPC2 = intersect(itPC,TheseTrials{ST},'stable');
%     
%     itPC3 = 1:length(TrPC);
%     itPC3(~ismember(itPC3,TheseTrials{ST})) = [];
%     [tPC4,itPC4] = sort(TrPC(itPC3),'descend');
    
    
    figure;
    subplot(4,1,1)
    plot(tPC,'-k')
    xlim([0.5 size(tPC,1)+0.5])
    set(gca,'xtick',[])
    title('stim 6')
    
    subplot(4,1,2:4)
    imagesc(CorrVals{ST}(:,itPC))
    colormap('gray')
%     colorbar
    set(gca,'clim',[0 1])
    
    orient(gcf,'landscape')
    print(gcf,fullfile(savedir,'stim6_act_sortPC'),'-dpdf','-bestfit')
    
    
    
    
    % Cluster population activity
    Data = CorrVals{ST}';
    
    % Calulate linkages
    Y = pdist(Data);     % cells X time
    Z = linkage(Y,'ward');
    
    leafOrder = optimalleaforder(Z,Y);
    
    % Plot dendrogram
    figure;
%     set(gcf,'Position',fullscreen)
    subplot(1,5,4:5);      % ,'ColorThreshold',40 ,'reorder',leafOrder
    [hd,tvals,outperm] = dendrogram(Z,0,'reorder',leafOrder,'ColorThreshold',30,'Orientation','right','Labels',cellfun(@num2str, num2cell([1:size(Data,1)]),'UniformOutput',false));
    set(gca,'tickdir','out','ytick',[])
    
    RespCluLabels = cluster(Z,'maxclust',3);
    
    
    % Plot trial activity
    subplot(1,5,2:3);
    
    imagesc( Data(fliplr(outperm),:) )
    colormap('gray')
%     colorbar
    set(gca,'clim',[0 1],'ytick',[])
%     set(gca,'ytick',1:size(Data,1),'yticklabel',fliplr(outperm))
    title('stim 6')
    xlabel('cells')
    
    % Get % correct this trial
    subplot(1,5,1);
    plot(TrPC(TheseTrials{ST}(outperm)),TheseTrials{ST},'-k')
    set(gca,'ytick',1:size(Data,1),'yticklabel',outperm)
    ylim([0.5 size(Data,1)+0.5])
    ylabel('trials')
    
    
    orient(gcf,'landscape')
    print(gcf,fullfile(savedir,'stim6_activity'),'-dpdf','-bestfit')
    
    
end

keyboard


%% Noise Correlation type stuff

% figure;
% xlim([-1 1])
% ylim([-1 1])
% hold on
% xlabel('Signal correlation')
% ylabel('Trial fidelity correlation')

for is = 1:numel(Sessions)
    
    % Get Subpop data
    datadir = fullfile(fn.figs,'ClassAM','AC','Full',Sessions{is});
    tabsavename = sprintf('CR_vFull_%s.mat',Sessions{is});
    q=load(fullfile(datadir,tabsavename));
    if size(q.CR,1)>1
        keyboard
    end
    
    ConfMat = mean(q.CR.Results{:},3);
    dpSP    = dp_from_ConfMat(ConfMat,crAmt);
    
    
    % Get SU data
    UnitIdx = q.CR.SUids{:}; %entries of CTTS/theseCells
    
    dpSU = nan(numel(UnitIdx),size(CTTS,4));
    for iu = 1:numel(UnitIdx)
        ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
        dpSU(iu,:) = dp_from_ConfMat(ConfMat,crAmt);
    end
    
    
    checkNoiseCorr(CTTS(UnitIdx,:,:,:));
    
    title(Sessions{is})
    text(-0.9,0.9,sprintf('%0.3f vs %0.3f\npop vs best SU', median(dpSP), median(max(dpSU,[],1)) ))
    
end

keyboard


% %     
% %     PSTHSize = TrainSize-TestSize;
% %     
% %     CorrVals = checkCorrAct(CTTS(UnitIdx,:,:,:),1,size(CTTS,4),numel(UnitIdx),'sim',PSTHSize);
% %     
% %     
% %     dpSU = nan(numel(UnitIdx),size(CTTS,4));
% %     for iu = 1:numel(UnitIdx)
% %         ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
% %         dpSU(iu,:) = dp_from_ConfMat(ConfMat,crAmt);
% %     end
% %     
% %     keyboard
% %     
% %     dpSP
% %     dpSU
% %     
% %     for ist = 1:size(CTTS,4)
% %         
% %         [rho,p] = corr(CorrVals{ist}');
% %         sum(sum(p<0.05))/2
% %         
% %         CorrVals{ist}
% %         
% %     end
% %     
    % Filter CTTS data
    
    
    
    
    % Plot tuning curves
%     hf2(is)=figure;    
%     plot(dpSU','LineWidth',2)
%     hold on
%     plot(dpSP,'k','LineWidth',4)
%     
%     xlim([0.8 8.2])    
%     xlabel('Stimuli')
%     ylabel('d''')
%     
%     title(sprintf('%s',Sessions{is}))
%     
%     print_eps_kp(hf2(is),fullfile(savedir,sprintf('dp_tuning_%s',Sessions{is})));
    
% end

keyboard
