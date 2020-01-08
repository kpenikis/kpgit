
close all

crAmt = 0.0001;

% Data settings
fn = set_paths_directories('','',1);
savedir = fullfile(fn.figs,'ClassAM','AC','Full','SessSubpops');

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

% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

StimOrder = [8 3 2 5 6 7 4 1];
newcolors = cmocean('thermal',length(StimOrder));


%% Total for stimulus set

% % figure;
% % hold on
% % 
% % for is = 1:numel(Sessions)
% %     
% %     datadir = fullfile(fn.figs,'ClassAM','AC','Full',Sessions{is});
% %     tabsavename = sprintf('CR_vFull_%s.mat',Sessions{is});
% %     
% %     q=load(fullfile(datadir,tabsavename));
% %     CR=q.CR;
% %     
% %     plot(CR.SUdps{:}',repmat(CR.dprime,[1 CR.iC]),...
% %         '+-','LineWidth',2,'Color',[0.1 0.1 0.1],...
% %         'MarkerSize',20)
% %     
% % end
% % 
% % plot([-0.2 2],[-0.2 2],':k')
% % xlim([-0.2 2])
% % ylim([-0.2 2])
% % axis square
% % 
% % xlabel('d'' of individual SU')
% % ylabel('d'' of subpopulation')
% % 
% % title('7 sessions, 10 mediocre SUs')
% % 
% % print_eps_kp(gcf,fullfile(savedir,'dp_SU_vs_Subpop'));
% % 
% % 
% % 
% % %% Each stimulus separately
% % 
% % ymax = 3.5;
% % 
% % for ist = StimOrder
% %     
% %     hf(ist)=figure;
% %     hold on
% %     
% %     for is = 1:numel(Sessions)
% %         
% %         % Get Subpop data
% %         datadir = fullfile(fn.figs,'ClassAM','AC','Full',Sessions{is});
% %         tabsavename = sprintf('CR_vFull_%s.mat',Sessions{is});
% %         q=load(fullfile(datadir,tabsavename));
% %         if size(q.CR,1)>1
% %             keyboard
% %         end
% %         
% %         ConfMat = mean(q.CR.Results{:},3);
% %         dpSP    = dp_from_ConfMat(ConfMat,crAmt,ist);
% %         
% %         UnitIdx = q.CR.SUids{:};
% %         
% %         % Get SU data
% %         
% %         dpSU = nan(numel(UnitIdx),1);
% %         for iu = 1:numel(UnitIdx)
% %             ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
% %             dpSU(iu) = dp_from_ConfMat(ConfMat,crAmt,ist);
% %         end
% %         
% %         % Plot it
% %         plot(dpSU',repmat(dpSP,[1 numel(UnitIdx)]),...
% %             '+-','LineWidth',2,'Color',newcolors(ist==StimOrder,:),...
% %             'MarkerSize',20)
% %         
% %     end
% %     
% %     plot([-0.2 ymax],[-0.2 ymax],':k')
% %     axis square
% %     
% %     xlim([-0.2 ymax])
% %     ylim([-0.2 ymax])
% %     
% %     xlabel('d'' of individual SU')
% %     ylabel('d'' of subpopulation')
% %     
% %     title(sprintf('stim %i',ist))
% %     
% %     print_eps_kp(hf(ist),fullfile(savedir,sprintf('dp_SU_vs_Subpop_st%i',ist)));
% %     
% % end
% % 
% % 
% % 
% % %% Within session d' tuning curves
% % 
% % for is = 1:numel(Sessions)
% %     
% %     % Get Subpop data
% %     datadir = fullfile(fn.figs,'ClassAM','AC','Full',Sessions{is});
% %     tabsavename = sprintf('CR_vFull_%s.mat',Sessions{is});
% %     q=load(fullfile(datadir,tabsavename));
% %     if size(q.CR,1)>1
% %         keyboard
% %     end
% %     
% %     ConfMat = mean(q.CR.Results{:},3);
% %     dpSP    = dp_from_ConfMat(ConfMat,crAmt);
% %     
% %     
% %     % Get SU data
% %     UnitIdx = q.CR.SUids{:};
% %     
% %     dpSU = nan(numel(UnitIdx),length(StimOrder));
% %     for iu = 1:numel(UnitIdx)
% %         ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
% %         dpSU(iu,:) = dp_from_ConfMat(ConfMat,crAmt);
% %     end
% %     
% %     % Plot it
% %     hf2(is)=figure;    
% %     plot(dpSU','LineWidth',2)
% %     hold on
% %     plot(dpSP,'k','LineWidth',4)
% %     
% %     xlim([0.8 8.2])    
% %     xlabel('Stimuli')
% %     ylabel('d''')
% %     
% %     title(sprintf('%s',Sessions{is}))
% %     
% %     print_eps_kp(hf2(is),fullfile(savedir,sprintf('dp_tuning_%s',Sessions{is})));
% %     
% % end



%% Max SU vs pop d' 

hf3=figure; 
hold on
plot([0.5 8.5],[0 0],':k')

threshDP=1;
nStClassified = nan(numel(Sessions),2);

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
    UnitIdx = q.CR.SUids{:};
    
    dpSU = nan(numel(UnitIdx),length(StimOrder));
    for iu = 1:numel(UnitIdx)
        ConfMat  = mean(CR_SU(UnitIdx(iu),:).Results{:},3);
        dpSU(iu,:) = dp_from_ConfMat(ConfMat,crAmt);
    end
    
    % Plot f3
    plot(1:8,dpSP-max(dpSU,[],1),'.','Color',[1 1 1]*0.12,'MarkerSize',20)
    
    % Save nStim > threshold 
    nStClassified(is,1) = sum(max(dpSU,[],1)>threshDP);
    nStClassified(is,2) = sum(dpSP>threshDP);
    
end

xlim([0.5 8.5])
xlabel('Stimulus')
ylabel('d'' Pop - best SU')


% print_eps_kp(hf3,fullfile(savedir,'BestSU_eachStim'));



% Now compare N stim classified
hf4=figure;
plot([0 8],[0 8],':k')
hold on

plot(nStClassified(:,1),nStClassified(:,2),'.','Color',[1 1 1]*0.10,'MarkerSize',30)

axis square
grid on
box off

set(gca,'xtick',1:8,'ytick',1:8,'Color','none')
xlim([0 8])
ylim([0 8])
xlabel('SUs')
ylabel('Subpop')
title(['N stimuli d''>' num2str(threshDP)])

print_eps_kp(hf4,fullfile(savedir,'NStim_thresh1'));


keyboard

