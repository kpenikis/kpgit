function CR_CatSum_Sess
% CR_CatSum_Sess
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
whichStim    = 'Speech';

crAmt = 0.01;


%%
% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,'CatSum','Sess');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

switch whichStim
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'AC'
        datadir = fullfile(fn.figs,'ClassAM','AC',whichClass,'Sess');
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'Speech'
        datadir = fullfile(fn.figs,'ClassSpeech','Speech',whichClass,'Sess');
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
end

% Load ALL CELLS results
% q=load(fullfile(datadir,'all',sprintf('CR_%s_%s.mat',whichClass,'all')));
% CR_Cat = q.CR;
% clear q
% q=load(fullfile(datadir,'all','Sum',sprintf('CR_%s_%s.mat',whichClass,'all')));
% CR_Sum = q.CR;
% clear q

% Load RS results
q=load(fullfile(datadir,'allRS',sprintf('CR_%s_%s.mat',whichClass,'allRS')));
CR_RS_Cat = q.CR;
clear q
q=load(fullfile(datadir,'allRS','Sum',sprintf('CR_%s_%s.mat',whichClass,'allRS')));
CR_RS_Sum = q.CR;
if strcmp(whichStim,'AC')
    CR_RS_Sum = CR_RS_Sum(2:end,:);
end
clear q

% Load NS results
% q=load(fullfile(datadir,'allNS',sprintf('CR_%s_%s.mat',whichClass,'allNS')));
% CR_NS_Cat = q.CR;
% clear q
% q=load(fullfile(datadir,'allNS','Sum',sprintf('CR_%s_%s.mat',whichClass,'allNS')));
% CR_NS_Sum = q.CR;
% clear q


% Recall which cells with sufficient trials
[~,theseCells] = recallDataParams(whichStim,'each',12);



%% Fig settings
set(groot,'DefaultTextInterpreter','tex')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
smallfig    = [1 scrsz(4)/3 scrsz(3)/4 scrsz(4)/3];
halfscreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

ymaxval =  5;

%==========================================================================
% RS indices
% iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
% iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);


%% Get individual stimulus d'

nStim    = size(CR_RS_Cat(1,:).Results{:},1);


% - - - All - - - 
% Cat
% dp_Cat = nan(size(CR_Cat,1),nStim);
% for ii=1:size(CR_Cat,1)    
%     dp_Cat(ii,:) = dp_from_ConfMat(mean(CR_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
% end
% [~,CI_Cat] = bootstrap4significance(CR_Cat);
% 
% % Sum
% dp_Sum = nan(size(CR_Sum,1),nStim);
% for ii=1:size(CR_Sum,1)    
%     dp_Sum(ii,:) = dp_from_ConfMat(mean(CR_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
% end
% [~,CI_Sum] = bootstrap4significance(CR_Sum);


% - - - RS - - - 
i_sim_RS = strcmp(CR_RS_Cat.trials,'sim');
% Cat
dp_RS_Cat = nan(size(CR_RS_Cat,1),nStim);
for ii=1:size(CR_RS_Cat,1)    
    dp_RS_Cat(ii,:) = dp_from_ConfMat(mean(CR_RS_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_RS_Cat] = bootstrap4significance(CR_RS_Cat);

% Sum
dp_RS_Sum = nan(size(CR_RS_Sum,1),nStim);
for ii=1:size(CR_RS_Sum,1)    
    dp_RS_Sum(ii,:) = dp_from_ConfMat(mean(CR_RS_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
end
[~,CI_RS_Sum] = bootstrap4significance(CR_RS_Sum);

% - - - NS - - - 
% Cat
% dp_NS_Cat = nan(size(CR_NS_Cat,1),nStim);
% for ii=1:size(CR_NS_Cat,1)    
%     dp_NS_Cat(ii,:) = dp_from_ConfMat(mean(CR_NS_Cat(ii,:).Results{:},3,'omitnan'),crAmt);
% end
% [~,CI_NS_Cat] = bootstrap4significance(CR_NS_Cat);
% 
% % Sum
% dp_NS_Sum = nan(size(CR_NS_Sum,1),nStim);
% for ii=1:size(CR_NS_Sum,1)    
%     dp_NS_Sum(ii,:) = dp_from_ConfMat(mean(CR_NS_Sum(ii,:).Results{:},3,'omitnan'),crAmt);
% end
% [~,CI_NS_Sum] = bootstrap4significance(CR_NS_Sum);





%%%%% PLOTS

%% Cat vs Sum for each stim, datapoint for each session 

figure;
set(gcf,'Position',smallfig)

plot([-2 6],[-2 6],'-k')
hold on
scatter(CR_RS_Cat(i_sim_RS,:).dprime,CR_RS_Sum(i_sim_RS,:).dprime,CR_RS_Sum(i_sim_RS,:).nC*6,'r','filled')
scatter(CR_RS_Cat(~i_sim_RS,:).dprime,CR_RS_Sum(~i_sim_RS,:).dprime,CR_RS_Sum(~i_sim_RS,:).nC*6,[0.5 0 0])
xlabel('Cat')
ylabel('Sum')
set(gca,'Color','none','xlim',[-0.5 3],'ylim',[-0.5 3])
axis square
title('Avg task')

print_eps_kp(gcf,fullfile(savedir,['CatSum_RS_' whichStim]))


hfo=figure;
set(hfo,'Position',halfscreen)
for ist = 1:nStim
    subplot(2,4,ist)
    hold on
    
    plot([-2 6],[-2 6],'-k')
    scatter(dp_RS_Cat(i_sim_RS,ist),dp_RS_Sum(i_sim_RS,ist),CR_RS_Sum(i_sim_RS,:).nC*6,'r','filled')
    scatter(dp_RS_Cat(~i_sim_RS,ist),dp_RS_Sum(~i_sim_RS,ist),CR_RS_Sum(~i_sim_RS,:).nC*6,[0.5 0 0])
    
    xlabel('Cat')
    ylabel('Sum')
    set(gca,'Color','none','xlim',[-2 6],'ylim',[-2 6])    
    axis square
    title(ist)
end
print_eps_kp(gcf,fullfile(savedir,['CatSum_RS_' whichStim '_each']))


%% Rand vs sim

figure;
set(gcf,'Position',smallfig)

plot([-2 6],[-2 6],'-k')
hold on

scatter(CR_RS_Sum(~i_sim_RS,:).dprime,CR_RS_Sum(i_sim_RS,:).dprime,CR_RS_Sum(i_sim_RS,:).nC*6,[50 205 50]./255,'filled')
scatter(CR_RS_Cat(~i_sim_RS,:).dprime,CR_RS_Cat(i_sim_RS,:).dprime,CR_RS_Cat(i_sim_RS,:).nC*6,[.14 .45 .49])
xlabel('Shuffled')
ylabel('Simultaneous')
set(gca,'Color','none','xlim',[-0.5 3],'ylim',[-0.5 3])
axis square
title('Avg task')

print_eps_kp(gcf,fullfile(savedir,['SimRand_RS_' whichStim]))


hfcat=figure;
set(hfcat,'Position',halfscreen)

for ist = 1:nStim
    subplot(2,4,ist)
    hold on
    
    plot([-2 6],[-2 6],'-k')
    scatter(dp_RS_Sum(~i_sim_RS,ist),dp_RS_Sum(i_sim_RS,ist),CR_RS_Sum(i_sim_RS,:).nC*6,[50 205 50]./255,'filled')
    scatter(dp_RS_Cat(~i_sim_RS,ist),dp_RS_Cat(i_sim_RS,ist),CR_RS_Cat(i_sim_RS,:).nC*6,[.14 .45 .49])
    
    xlabel('Shuffled')
    ylabel('Simultaneous')
    set(gca,'Color','none','xlim',[-2 6],'ylim',[-2 6])    
    axis square
    title(ist)
end
print_eps_kp(gcf,fullfile(savedir,['SimRand_RS_' whichStim '_each']))



%%

keyboard

hf=figure;
set(hf,'Position',halfscreen)

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


% Difference for each stimulus
subplot(2,2,2)
hold on
bar(0+(1:numel(dp_Cat)),dp_Cat-dp_Sum,1,'k')
bar(9+(1:numel(dp_RS_Cat)),dp_RS_Cat-dp_RS_Sum,1,'r')
bar(18+(1:numel(dp_NS_Cat)),dp_NS_Cat-dp_NS_Sum,1,'b')

set(gca,'Color','none','xlim',[0 27],'ylim',[-1.5 2.5],...
    'xtick',[1:numel(dp_RS_Cat) 9+(1:numel(dp_NS_Cat)) 18+(1:numel(dp_NS_Cat))],...
    'xticklabel',[1:8 1:8 1:8])
ylabel('\Delta d''')
xlabel('Stimulus')

subplot(2,2,4)
hold on
plot(1:numel(dp_Cat),dp_Cat,'k--','LineWidth',2)
plot(1:numel(dp_Sum),dp_Sum,'k-','LineWidth',2)
plot(9+(1:numel(dp_RS_Cat)),dp_RS_Cat,'r--','LineWidth',2)
plot(9+(1:numel(dp_RS_Sum)),dp_RS_Sum,'r-','LineWidth',2)
plot(18+(1:numel(dp_NS_Cat)),dp_NS_Cat,'b--','LineWidth',2)
plot(18+(1:numel(dp_NS_Cat)),dp_NS_Sum,'b-','LineWidth',2)

set(gca,'Color','none','xlim',[0 27],'ylim',[0 ymaxval],...
    'xtick',[1:numel(dp_RS_Cat) 9+(1:numel(dp_NS_Cat)) 18+(1:numel(dp_NS_Cat))],...
    'xticklabel',[1:8 1:8 1:8])
ylabel('d''')
xlabel('Stimulus')



print_eps_kp(gcf,fullfile(savedir,['AllRSNS_' whichStim]))


end








