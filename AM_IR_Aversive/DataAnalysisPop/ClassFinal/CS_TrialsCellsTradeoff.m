% compare ActVec and Projections classifiers

% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',12)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');
% New Results, with CReach run first

scrsz = get(0,'ScreenSize');     %[left bottom width height]
widescrn  = [1 scrsz(4)/3 scrsz(3) scrsz(4)/3];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];

close all
whichStim = 'Speech';
whichClass = 'Full';


dpm = 7;

fn = set_paths_directories;
savedir = fullfile(fn.figs,'ClassResults','Corrections','TrsCellsTradeoff');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
% results for adding RS cells dp ranked

fn = set_paths_directories;
switch whichStim
    case 'AC'
        rootdir = fullfile(fn.figs,'ClassAM',whichStim);
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassSpeech',whichStim);
end
Ntrs = dir(fullfile(rootdir,whichClass,'minTrs*'));

Res=struct;
for ii = 1:numel(Ntrs)
    
    clear CR
    load(fullfile(rootdir,whichClass,Ntrs(ii).name,'dpRank_RS',['CR_v' whichClass '_dpRank_RS.mat']))
    
    Res(ii).minTrs   = str2num(Ntrs(ii).name(end-1:end));
    [Res(ii).NC,inc] = sort(CR.nC);
    Res(ii).dp       = CR.dprime(inc);
    Res(ii).PC       = CR.PC(inc);
    
    % Also load CReach results
    switch whichClass
        case {'Full' 'Sum'}
            clear CR
            load(fullfile(rootdir,'Full',Ntrs(ii).name,'CR_each.mat'))
            Res(ii).bestSU   = max(CR.dprime);
            Res(ii).bestSUpc = max(CR.PC);
        case {'ActVec' 'SumAV'}
            clear CR
            load(fullfile(rootdir,'ActVec',Ntrs(ii).name,'CR_each.mat'))
            Res(ii).bestSU   = max(CR.dprime);
            Res(ii).bestSUpc = max(CR.PC);
    end
    
end


%% Raw results (d' as a fct of N cells, for each trial minimum cutoff

hf=figure;
hold on
set(gcf,'Position',widescrn)

xlabelpos=[];
xlabelstr={};
clear ip
for imt = 1:numel(Res)
    
    ip(imt) = plot([1; Res(imt).NC]+200*(imt-1),min([Res(imt).bestSU;Res(imt).dp],dpm),'LineWidth',3);
    
%     xlabelpos = [xlabelpos; 0+200*(imt-1); Res(imt).NC+200*(imt-1)];
%     xlabelstr = [xlabelstr; sprintf('%i Trs',Res(imt).minTrs); cellstr(num2str(Res(imt).NC)) ];
    xlabelpos = [xlabelpos; max(Res(imt).NC)+200*(imt-1)];
    xlabelstr = [xlabelstr; cellstr(num2str(max(Res(imt).NC))) ];    
end

legend(ip,cellfun(@num2str, {Res.minTrs},'UniformOutput',0),'Location','best')

xlabel('Total # Cells in ensemble')
ylabel('d''')
set(gca,'Color','none','xtick',xlabelpos,'xticklabel',xlabelstr)
ylim([0 dpm])

savename = sprintf('%s_%s_dpAddCells_dpm%i',whichStim,whichClass,dpm);
print_eps_kp(hf,fullfile(savedir,savename))



%% Visualize trade off

% collect d's
dp_allRS = nan(numel(Res),1);
dp_max   = nan(numel(Res),1);
dp_SU    = nan(numel(Res),1);
Ncells   = nan(numel(Res),1);
for imt = 1:numel(Res)
    dp_allRS(imt) = min(Res(imt).dp(end),dpm);
    dp_max(imt)   = min(max(Res(imt).dp),dpm);
    dp_SU(imt)    = min(max(Res(imt).bestSU),dpm);
    Ncells(imt)   = Res(imt).NC(end);
end

% sort by d' with all cells
% [~,inmt] = sort(dp_allRS);
inmt = [1:numel(Res)]';

% plot min N trs by [d'all - d'30]
% label d' max
hf2 = figure;
set(gcf,'Position',widesmall)

subplot(1,2,1)
% plot([0 500],[0 0],'--k')
hold on
clear ip2
for iinmt = inmt'
%     plot(Res(iinmt).minTrs, dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50)
    ip2(iinmt==inmt)=plot(Ncells(iinmt), dp_allRS(iinmt)-dp_max(iinmt),'.','MarkerSize',50);
end

legend(ip2,cellfun(@num2str, {Res(inmt).minTrs},'UniformOutput',0),'Location','best')
set(gca,'Color','none')
xlabel('maximum N cells')
xlim([0 200])
ylabel('d''_all  -  max d''')


subplot(1,2,2)
% plot([0 500],[0 0],'--k')
hold on
clear ip2
for iinmt = inmt'
%     plot(Res(iinmt).minTrs, dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50)
    ip2(iinmt==inmt)=plot(Ncells(iinmt), dp_allRS(iinmt)-dp_SU(iinmt),'.','MarkerSize',50);
end

legend(ip2,cellfun(@num2str, {Res(inmt).minTrs},'UniformOutput',0),'Location','best')
set(gca,'Color','none')
xlabel('maximum N cells')
xlim([0 200])
ylabel('d''_all  -  best SU d''')


savename = sprintf('%s_%s_NcellsDPdiff',whichStim,whichClass);
print_eps_kp(hf2,fullfile(savedir,savename))


%%  PERCENT CORRECT


%% Raw results (d' as a fct of N cells, for each trial minimum cutoff

hf3=figure;
hold on
set(gcf,'Position',widescrn)

xlabelpos=[];
xlabelstr={};
clear ip
for imt = 1:numel(Res)
    
    ip(imt) = plot([1; Res(imt).NC]+200*(imt-1),[Res(imt).bestSUpc;Res(imt).PC],'LineWidth',3);
    
%     xlabelpos = [xlabelpos; 0+200*(imt-1); Res(imt).NC+200*(imt-1)];
%     xlabelstr = [xlabelstr; sprintf('%i Trs',Res(imt).minTrs); cellstr(num2str(Res(imt).NC)) ];
    xlabelpos = [xlabelpos; max(Res(imt).NC)+200*(imt-1)];
    xlabelstr = [xlabelstr; cellstr(num2str(max(Res(imt).NC))) ];    
end

legend(ip,cellfun(@num2str, {Res.minTrs},'UniformOutput',0),'Location','best')

xlabel('Total # Cells in ensemble')
ylabel('Percent Correct')
set(gca,'Color','none','xtick',xlabelpos,'xticklabel',xlabelstr)
ylim([0 100])

savename = sprintf('%s_%s_dpAddCells_PC',whichStim,whichClass);
print_eps_kp(hf3,fullfile(savedir,savename))



%% Visualize trade off

% collect PCs
PC_allRS = nan(numel(Res),1);
PC_max   = nan(numel(Res),1);
PC_SU    = nan(numel(Res),1);
Ncells   = nan(numel(Res),1);
for imt = 1:numel(Res)
    PC_allRS(imt) = Res(imt).PC(end);
    PC_max(imt)   = max(Res(imt).PC);
    PC_SU(imt)    = max(Res(imt).bestSUpc);
    Ncells(imt)   = Res(imt).NC(end);
end

% sort by d' with all cells
% [~,inmt] = sort(dp_allRS);
inmt = [1:numel(Res)]';

% plot min N trs by [d'all - d'30]
% label d' max
hf4 = figure;
set(gcf,'Position',widesmall)

subplot(1,2,1)
% plot([0 500],[0 0],'--k')
hold on
clear ip2
for iinmt = inmt'
%     plot(Res(iinmt).minTrs, dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50)
    ip2(iinmt==inmt)=plot(Ncells(iinmt), PC_allRS(iinmt)-PC_max(iinmt),'.','MarkerSize',50);
end

legend(ip2,cellfun(@num2str, {Res(inmt).minTrs},'UniformOutput',0),'Location','best')
set(gca,'Color','none')
xlabel('maximum N cells')
xlim([0 200])
ylabel('PC_all  -  max PC')


subplot(1,2,2)
% plot([0 500],[0 0],'--k')
hold on
clear ip2
for iinmt = inmt'
%     plot(Res(iinmt).minTrs, dp_allRS(iinmt)-dp_30RS(iinmt),'.','MarkerSize',50)
    ip2(iinmt==inmt)=plot(Ncells(iinmt), PC_allRS(iinmt)-PC_SU(iinmt),'.','MarkerSize',50);
end

legend(ip2,cellfun(@num2str, {Res(inmt).minTrs},'UniformOutput',0),'Location','best')
set(gca,'Color','none')
xlabel('maximum N cells')
xlim([0 200])
ylabel('PC_all  -  best SU PC')


savename = sprintf('%s_%s_NcellsPCdiff',whichStim,whichClass);
print_eps_kp(hf4,fullfile(savedir,savename))








