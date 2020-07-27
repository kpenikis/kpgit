% compare ActVec and Projections classifiers

% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',12)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

%%
% results for adding RS cells dp ranked

close all

fn = set_paths_directories;

% Load original: min 12 trials
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/allRS/CR_vFull_allRS.mat')
CR_all_Proj = CR;
clear CR
CR_all_Proj = CR_all_Proj(1,:);

load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/AC/Full/dpRank_RS/CR_vFull_dpRank_RS.mat')
CR_Proj = CR;
clear CR

CR_Proj = CR_Proj(CR_Proj.iC==1,:);
[nC_P,inC_P] = sort(CR_Proj.nC);
CR_Proj = CR_Proj(inC_P,:);

Res=struct;
Res(1).minTrs = 12;
Res(1).NC     = [nC_P; 181];
Res(1).dp     = [CR_Proj.dprime; CR_all_Proj.dprime];


% Load results with other min N trials
rootdir = fullfile(fn.figs,'ClassAM','AC','Full','dpRank_RS');

Ntrials = [22 23 24 27 32];

for imt = 1:numel(Ntrials)
    
    
    load(fullfile(rootdir,['Trs' num2str(Ntrials(imt))],'CR_vFull_dpRank_RS.mat'))
    
    Res(imt+1).minTrs = Ntrials(imt);
    Res(imt+1).NC     = CR.nC;
    Res(imt+1).dp     = CR.dprime;
    
    
end



%% Raw results (d' as a fct of N cells, for each trial minimum cutoff

hf=figure;
hold on

xlabelpos=[];
xlabelstr={};
for imt = 1:numel(Res)
    
    ip(imt) = plot(Res(imt).NC+200*(imt-1),Res(imt).dp,'LineWidth',2);
    
    xlabelpos = [xlabelpos; 0+200*(imt-1); Res(imt).NC+200*(imt-1)];
    xlabelstr = [xlabelstr; sprintf('%i Trs',Res(imt).minTrs); cellstr(num2str(Res(imt).NC)) ];
    
end

legend(ip,cellfun(@num2str, {Res.minTrs},'UniformOutput',0),'Location','best')

xlabel('# Cells in ensemble')
ylabel('d''')
set(gca,'Color','none','xtick',xlabelpos,'xticklabel',xlabelstr)
ylim([0 5])


savedir = fullfile(fn.figs,'ClassResults','Corrections');
print_eps_kp(hf,fullfile(savedir,'TCtradeoff_addCells'))


%% Visualize trade off

% collect d's
dp_allRS = nan(numel(Res),1);
dp_30RS  = nan(numel(Res),1);
Ncells   = nan(numel(Res),1);
for imt = 1:numel(Res)
    dp_allRS(imt) = Res(imt).dp(end);
    dp_30RS(imt)  = Res(imt).dp(Res(imt).NC==30);
    Ncells(imt)   = Res(imt).NC(end);
end

% sort by d' with all cells
[~,inmt] = sort(dp_allRS);

% plot min N trs by [d'all - d'30]
% label d' max
hf2 = figure;
plot([0 500],[0 0],'--k')
hold on
for iinmt = inmt'
%     plot(Res(iinmt).minTrs, dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50)
    ip2(iinmt==inmt)=plot(Ncells(iinmt), dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50);
end

legend(ip2,cellfun(@num2str, {Res(inmt).minTrs},'UniformOutput',0),'Location','best')
set(gca,'Color','none')
xlabel('maximum N cells')
xlim([0 200])
ylabel('d''_all  -  d''_30cells')


print_eps_kp(hf2,fullfile(savedir,'TCtradeoff_NcellsDPdiff'))


