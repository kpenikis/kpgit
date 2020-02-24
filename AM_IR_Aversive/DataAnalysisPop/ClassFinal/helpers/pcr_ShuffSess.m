

whichStim = 'Speech';
whichCells = 'RS';

% Data settings
fn = set_paths_directories('','',1);

switch whichStim
    case 'AC'    
        keyboard
        
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % ~~~~~~ SHAPE ~~~~~~
        % Full class results
        datadir = fullfile(fn.figs,'ClassAM',whichStim,'TrShuff','SessionsMin2cells',whichCells);
        q=load(fullfile(datadir,'CR_vTrShuff.mat'));
        CR_Sh_Full = q.CR;
        clear q
        
        % Rate class results
        datadir = fullfile(fn.figs,'ClassAM',whichStim,'TrShuff_Rate','SessionsMin2cells',whichCells);
        q=load(fullfile(datadir,'CR_vTrShuff_Rate.mat'));
        CR_Sh_Rate = q.CR;
        clear q
        
        % ~~~~~ CONTEXT ~~~~~
%         datadir = fullfile(fn.figs,'ClassContext','AC',whichClass,'each');
%         q=load(fullfile(datadir,'CR_each.mat'));
%         CRe_Context = q.CR;
%         clear q
        
        savedir = fullfile(fn.figs,'ClassResults','TrialShuffle',whichStim);
        
    case 'Speech'        
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % ~~~~~~ SHAPE ~~~~~~
        % Temp class results
        datadir = fullfile(fn.figs,'ClassSpeech',whichStim,'TrShuff','SessionsMin3cells');
        q=load(fullfile(datadir,'CR_vTrShuff.mat'));
        CR_Sh_Temp = q.CR;
        clear q
        
        % Rate class results
        datadir = fullfile(fn.figs,'ClassSpeech',whichStim,'TrShuff_Rate','SessionsMin2cells',whichCells);
        q=load(fullfile(datadir,'CR_vTrShuff_Rate.mat'));
        CR_Sh_Rate = q.CR;
        clear q
        
        % ~~~~~ CONTEXT ~~~~~
%         datadir = fullfile(fn.figs,'ClassContext','AC',whichClass,'each');
%         q=load(fullfile(datadir,'CR_each.mat'));
%         CRe_Context = q.CR;
%         clear q
        
        savedir = fullfile(fn.figs,'ClassResults','TrialShuffle',whichStim);
end

if ~exist(savedir,'dir')
    mkdir(savedir)
end


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];



%% 

CellMin = 3;

CR_Temp = CR_Sh_Temp(CR_Sh_Temp.iC>CellMin,:);
CR_Rate = CR_Sh_Rate(CR_Sh_Rate.iC>CellMin,:);


% Shuffled vs simultaneous trials

% Temporal code
figure;
plot([0 3],[0 3],'-k')
hold on
scatter([CR_Temp.dprime(1:2:end)],[CR_Temp.dprime(2:2:end)],20.*CR_Temp.iC(1:2:end),'k','filled')
xlabel('d'' shuffled trials')
ylabel('d'' simultaneous trials')
axis square
title('Temporal')

print_eps_kp(gcf,fullfile(savedir,['TrShuff_' whichCells '_Temp2']))


% Rate code
figure;
plot([-0.2 1],[-0.2 1],'-k')
hold on
scatter([CR_Rate.dprime(1:2:end)],[CR_Rate.dprime(2:2:end)],20.*CR_Rate.iC(1:2:end),'k','filled')
xlabel('d'' shuffled trials')
ylabel('d'' simultaneous trials')
axis square
title('Rate')
xlim([-0.2 1])
ylim([-0.2 1])

print_eps_kp(gcf,fullfile(savedir,['TrShuff_' whichCells '_Rate2']))




