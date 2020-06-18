% NOT yet broken down by stimulus 

close all

whichClass   = 'ActVec';

crAmt = 0.01;

% Data settings
fn = set_paths_directories('','',1);

savedir = fullfile(fn.figs,'ClassResults',whichClass);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~

datadir_AM = fullfile(fn.figs,'ClassAM','AC',whichClass);
q = load(fullfile(fn.processed,'Units'));
UnitInfo_AM = q.UnitInfo;
UnitData_AM = q.UnitData;
clear q

% Load SU data 
q=load(fullfile(datadir_AM,'each','CR_each.mat'));
CReach_AM = q.CR;
clear q

% Load subpop results    
% q=load(fullfile(datadir_AM,'pkFR_RS',['CR_v' whichClass '_pkFR_RS.mat']));
% CR_pfr_AM = q.CR;
% clear q
% q=load(fullfile(datadir_AM,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
% CR_Qpfr_AM = q.CR;
% clear q
q=load(fullfile(datadir_AM,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']));
CR_dp_AM = q.CR;
clear q
q=load(fullfile(datadir_AM,'maxdp_RS',['CR_v' whichClass '_maxdp_RS.mat']));
CR_mdp_AM = q.CR;
clear q


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
sqsmall   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];


%%

fstCells = unique(CR_mdp_AM.iC);
numCells = unique(CR_mdp_AM.nC);

FCcolors = cmocean('phase',4);

figure;
plot([0 5],[0 5],'k')
hold on
for iC = 1:numel(fstCells)
    for nC = 1:numel(numCells)
        imxdp = find(CR_mdp_AM.iC==fstCells(iC) & CR_mdp_AM.nC==numCells(nC));
        iavdp = find(CR_dp_AM.iC==fstCells(iC) & CR_dp_AM.nC==numCells(nC));
        
        if ~isempty(iavdp) && ~isempty(imxdp)
            plot(CR_dp_AM.dprime(iavdp),CR_mdp_AM.dprime(imxdp),...
                'o','MarkerSize',10,'MarkerFaceColor',FCcolors(iC,:),'MarkerEdgeColor','none')
        end
    end
end
axis square
xlabel('pool d'' from avg dp rank')
ylabel('pool d'' from MAX dp rank')

print_eps_kp(gcf,fullfile(savedir,'Compare_pools_byMaxbyAvg'))









