function CR_SUvPop_overtime
% CR_SUvPop_overtime
%
%   Compare classification when cells are kept separate or when activity
%   summed. RS and NS cells separately, and all cells together.
%
%   Classification using the entire 500 ms segment.
%
% KP, 2020-04
%

close all

whichClass   = 'Full';
whichStim    = 'AC';
WinFolder    = 'Sliding100';

crAmt = 0.01;


%%
% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,WinFolder);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

switch whichStim
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'AC'
        datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass);
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo_AM = q.UnitInfo;
        UnitData_AM = q.UnitData;
        clear q
        
        % Load SU data
        q=load(fullfile(datadir_AM,'each',WinFolder,'CR_each.mat'));
        CRe = q.CR;
        clear q
        
        
        % Load Pop results
        
        % All RS
        q=load(fullfile(datadir_AM,'allRS',WinFolder,sprintf('CR_v%s_%s.mat',whichClass,'allRS')));
        CR_RS = q.CR;
        clear q
        
        % All NS
        q=load(fullfile(datadir_AM,'allNS',WinFolder,sprintf('CR_v%s_%s.mat',whichClass,'allNS')));
        CR_NS = q.CR;
        clear q
        
        
        % Recall which cells had sufficient trials
        [~,theseCells_AM] = recallDataParams('AC','each',12);
        
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'Speech'
        % datadir_Sp = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
        % q = load(fullfile(fn.processed,'UnitsVS'));
        % UnitInfo_Sp = q.UnitInfo;
        % UnitData_Sp = q.UnitData;
        % clear q
        %
end


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

ymaxval =  3;


%% Get a few variables

% RS indices
iRS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak>0.43);
iNS_AM = find(UnitInfo_AM(theseCells_AM,:).TroughPeak<=0.43);
% iRS_Sp = find(UnitInfo_Sp(theseCells_Sp,:).TroughPeak>0.43);

% Time vector
[x_RS ,ixRS]= sort(CR_RS.WinEnd-550);
[x_NS ,ixNS]= sort(CR_NS.WinEnd-550);


%%

hf=figure;
set(hf,'Position',halfscreen)

subplot(1,2,1)

plot(x_RS,CR_RS.dprime(ixRS),'r-','LineWidth',2)
hold on
plot(x_NS,CR_NS.dprime(ixNS),'b-','LineWidth',2)

% Find SU data matching windows
x_SU = unique(CRe.WinEnd);
mdp  = nan(size(x_SU));
mndp = nan(size(x_SU));
for ii = 1:numel(x_SU)
    plotSpread(CRe(CRe.WinEnd==x_SU(ii),:).dprime,...
        'distributionIdx',x_SU(ii)-550*ones(sum(CRe.WinEnd==x_SU(ii)),1),...
        'distributionColor','k')
    
    mdp(ii)  = median(CRe(CRe.WinEnd==x_SU(ii),:).dprime);
    mndp(ii) = mean(CRe(CRe.WinEnd==x_SU(ii),:).dprime);
end
plot(x_SU-550,mdp,'Color',[1 0.8 0],'LineWidth',2)


set(gca,'Color','none','xlim',[0 500],'ylim',[0 ymaxval],...
    'xtick',0:100:500,'xticklabel',0:100:500)
xlabel('Middle of analysis window (ms)')
ylabel('d''')
title(WinFolder)

print_eps_kp(hf,fullfile(savedir,'allRS-allNS-SU'))


keyboard
end








