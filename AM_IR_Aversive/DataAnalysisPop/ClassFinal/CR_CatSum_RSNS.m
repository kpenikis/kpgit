function CR_CatSum_RSNS
% CR_CatSum_RSNS
%  fig 10
%
%   Compare classification when cells are kept separate or when activity
%   summed. RS and NS cells separately, and all cells together.
%
%   Classification using the entire 500 ms segment.
%
% KP, 2020-04
%

% close all

whichClass   = 'Full';
whichStim    = 'AC';

crAmt = 0.01;


%%
% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,'CatSum');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

switch whichStim
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'AC'
        datadir = fullfile(fn.figs,'ClassAM','AC',whichClass);
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'Speech'
        datadir = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
end

% Load ALL CELLS results
q=load(fullfile(datadir,'all',sprintf('CR_v%s_%s.mat',whichClass,'all')));
CR_Cat = q.CR;
clear q
q=load(fullfile(datadir,'all','Sum',sprintf('CR_v%s_%s.mat',whichClass,'all')));
CR_Sum = q.CR;
clear q

% Load RS results
q=load(fullfile(datadir,'allRS',sprintf('CR_v%s_%s.mat',whichClass,'allRS')));
CR_RS_Cat = q.CR;
clear q
q=load(fullfile(datadir,'allRS','Sum',sprintf('CR_v%s_%s.mat',whichClass,'allRS')));
CR_RS_Sum = q.CR;
clear q

% Load NS results
q=load(fullfile(datadir,'allNS',sprintf('CR_v%s_%s.mat',whichClass,'allNS')));
CR_NS_Cat = q.CR;
clear q
q=load(fullfile(datadir,'allNS','Sum',sprintf('CR_v%s_%s.mat',whichClass,'allNS')));
CR_NS_Sum = q.CR;
clear q


% Recall which cells with sufficient trials
[~,theseCells] = recallDataParams(whichStim,'each',12);

switch whichStim
    case 'AC'
        theseStim  = 1:8;
    case 'DB'
        theseStim  = [1:6 9:10];
    case 'Speech'
        theseStim  = [4 3 2 1 5 6 7 8]; %1:size(Cell_Time_Trial_Stim,4);
end


%% Fig settings
set(groot,'DefaultTextInterpreter','tex')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];
qrtscreen   = [1 scrsz(4)/4 scrsz(3) scrsz(4)/4];

ymaxval =  5;

%==========================================================================
% RS indices
% iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
% iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);


%% Get individual stimulus d'

nStim    = size(CR_RS_Cat(1,:).Results{:},1);


% - - - All - - - 
% Cat
dp_Cat = nan(size(CR_Cat,1),nStim);
for ii=1:size(CR_Cat,1)    
    dp_Cat(ii,theseStim) = dp_from_ConfMat(mean(CR_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_Cat] = bootstrap4significance(CR_Cat);

% Sum
dp_Sum = nan(size(CR_Sum,1),nStim);
for ii=1:size(CR_Sum,1)    
    dp_Sum(ii,theseStim) = dp_from_ConfMat(mean(CR_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_Sum] = bootstrap4significance(CR_Sum);


% - - - RS - - - 
% Cat
CR_RS_Cat = CR_RS_Cat(~CR_RS_Cat.exclSpec & ~CR_RS_Cat.exNonSig,:);
dp_RS_Cat = nan(size(CR_RS_Cat,1),nStim);
for ii=1:size(CR_RS_Cat,1)    
    dp_RS_Cat(ii,theseStim) = dp_from_ConfMat(mean(CR_RS_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_RS_Cat] = bootstrap4significance(CR_RS_Cat);

% Sum
dp_RS_Sum = nan(size(CR_RS_Sum,1),nStim);
for ii=1:size(CR_RS_Sum,1)    
    dp_RS_Sum(ii,theseStim) = dp_from_ConfMat(mean(CR_RS_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_RS_Sum] = bootstrap4significance(CR_RS_Sum);

% - - - NS - - - 
% Cat
dp_NS_Cat = nan(size(CR_NS_Cat,1),nStim);
for ii=1:size(CR_NS_Cat,1)    
    dp_NS_Cat(ii,theseStim) = dp_from_ConfMat(mean(CR_NS_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_NS_Cat] = bootstrap4significance(CR_NS_Cat);

% Sum
dp_NS_Sum = nan(size(CR_NS_Sum,1),nStim);
for ii=1:size(CR_NS_Sum,1)    
    dp_NS_Sum(ii,theseStim) = dp_from_ConfMat(mean(CR_NS_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_NS_Sum] = bootstrap4significance(CR_NS_Sum);



%% Just RS cells

hf=figure;
set(hf,'Position',qrtscreen)

% Task average
subplot(1,5,1)
hold on

%RS
bar(1,CR_RS_Cat.dprime,'r')
plot([4 4],CI_RS_Cat,'-k','LineWidth',2)

bar(2,CR_RS_Sum.dprime,'r')
plot([5 5],CI_RS_Sum,'-k','LineWidth',2)

set(gca,'Color','none','xlim',[0 3],'ylim',[0 ymaxval],...
    'xtick',[1 2],'xticklabel',{'RS cat' 'RS sum'})
ylabel('d''')


%%% Set order of stimuli
[~,iStimRank] = sort(dp_RS_Cat-dp_RS_Sum);

% Difference for each stimulus
subplot(1,5,4:5)
bar(0+(1:numel(dp_RS_Cat)),dp_RS_Cat(iStimRank)-dp_RS_Sum(iStimRank),1,'r')

set(gca,'Color','none','xlim',[0 9],'ylim',[-1.5 2.5],...
'xtick',1:numel(dp_RS_Cat),'xticklabel',iStimRank)
ylabel('\Delta d''')
xlabel('Stimulus')

% Line for each stimulus
subplot(1,5,2:3)
hold on
plot(0+(1:numel(dp_RS_Cat)),dp_RS_Cat(iStimRank),'r--','LineWidth',2)
plot(0+(1:numel(dp_RS_Sum)),dp_RS_Sum(iStimRank),'r-','LineWidth',2)

set(gca,'Color','none','xlim',[0 9],'ylim',[0 ymaxval],...
    'xtick',1:numel(dp_RS_Cat),'xticklabel',iStimRank)
ylabel('d''')
xlabel('Stimulus')



print_eps_kp(gcf,fullfile(savedir,['RS_ranked_' whichStim]))



%% All, RS, & NS

keyboard

hf=figure;
set(hf,'Position',halfscreen)

% Task average
subplot(2,2,1)
hold on

bar(1,CR_Cat.dprime,'k')
plot([1 1],CI_Cat,'-k','LineWidth',2)

bar(2,CR_Sum.dprime,'k')
plot([2 2],CI_Sum,'-k','LineWidth',2)

%RS
bar(4,CR_RS_Cat.dprime,'r')
plot([4 4],CI_RS_Cat,'-k','LineWidth',2)

bar(5,CR_RS_Sum.dprime,'r')
plot([5 5],CI_RS_Sum,'-k','LineWidth',2)

%NS
bar(7,CR_NS_Cat.dprime,'b')
plot([7 7],CI_NS_Cat,'-k','LineWidth',2)

bar(8,CR_NS_Sum.dprime,'b')
plot([8 8],CI_NS_Sum,'-k','LineWidth',2)

set(gca,'Color','none','xlim',[0 9],'ylim',[0 ymaxval],...
    'xtick',[1 2 4 5 7 8],'xticklabel',{'All cat' 'All sum' 'RS cat' 'RS sum' 'NS cat' 'NS sum'})
ylabel('d''')


%%% Set order of stimuli
[~,iStimRank] = sort(dp_RS_Cat-dp_RS_Sum);

% Difference for each stimulus
subplot(2,2,4)
hold on
bar(0+(1:numel(dp_Cat)),dp_Cat(iStimRank)-dp_Sum(iStimRank),1,'k')
bar(9+(1:numel(dp_RS_Cat)),dp_RS_Cat(iStimRank)-dp_RS_Sum(iStimRank),1,'r')
bar(18+(1:numel(dp_NS_Cat)),dp_NS_Cat(iStimRank)-dp_NS_Sum(iStimRank),1,'b')

set(gca,'Color','none','xlim',[0 27],'ylim',[-1.5 2.5],...
    'xtick',[1:numel(dp_RS_Cat) 9+(1:numel(dp_NS_Cat)) 18+(1:numel(dp_NS_Cat))],...
    'xticklabel',[(iStimRank) (iStimRank) (iStimRank)])
ylabel('\Delta d''')
xlabel('Stimulus')

% Line for each stimulus
subplot(2,2,2)
hold on
plot(1:numel(dp_Cat),dp_Cat(iStimRank),'k--','LineWidth',2)
plot(1:numel(dp_Sum),dp_Sum(iStimRank),'k-','LineWidth',2)
plot(9+(1:numel(dp_RS_Cat)),dp_RS_Cat(iStimRank),'r--','LineWidth',2)
plot(9+(1:numel(dp_RS_Sum)),dp_RS_Sum(iStimRank),'r-','LineWidth',2)
plot(18+(1:numel(dp_NS_Cat)),dp_NS_Cat(iStimRank),'b--','LineWidth',2)
plot(18+(1:numel(dp_NS_Cat)),dp_NS_Sum(iStimRank),'b-','LineWidth',2)

set(gca,'Color','none','xlim',[0 27],'ylim',[0 ymaxval],...
    'xtick',[1:numel(dp_RS_Cat) 9+(1:numel(dp_NS_Cat)) 18+(1:numel(dp_NS_Cat))],...
    'xticklabel',[(iStimRank) (iStimRank) (iStimRank)])
ylabel('d''')
xlabel('Stimulus')



print_eps_kp(gcf,fullfile(savedir,['AllRSNS_ranked_' whichStim]))




%% bar for each stim

% hfs=figure;
% set(hfs,'Position',thrdscreen)
% for ist = 1:nStim
%     subplot(2,4,ist)
%     hold on
%     bar(1,dp_RS_Cat(:,ist),'r')
%     
%     bar(2,dp_RS_Sum(:,ist),'r')
%     
%     bar(4,dp_NS_Cat(:,ist),'b')
%     
%     bar(5,dp_NS_Sum(:,ist),'b')
%     
%     set(gca,'Color','none','xlim',[0 6],'ylim',[0 ymaxval],...
%         'xtick',[1 2 4 5],'xticklabel',{'cat' 'sum' 'cat' 'sum'})    
%     ylabel('d''')
%     title(ist)
% end


end








