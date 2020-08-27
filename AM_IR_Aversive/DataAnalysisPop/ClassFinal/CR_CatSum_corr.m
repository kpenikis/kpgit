function CR_CatSum_corr
% Plots Fig 10* last panel
% new data

% close all

whichClass   = 'Full';
whichCells   = 'dpRank_RS';

whichStim    = 'Speech';

minNtrs      = 18;

yaxvalue     = 'dp';

dpmaxval     = 4;
crAmt        = 0.001;

%%

% Data settings
fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults','Corrections','CatSum','Thesis');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


switch whichStim
    %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'AC'
        
        rootdir = fullfile(fn.figs,'ClassAM',whichStim,whichClass,['minTrs' num2str(minNtrs)]);
        switch whichClass
            case 'Full'
                rootdir_sum = fullfile(fn.figs,'ClassAM',whichStim,'Sum',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'Sum';
            case 'ActVec'
                rootdir_sum = fullfile(fn.figs,'ClassAM',whichStim,'SumAV',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'SumAV';
        end
        
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        theseStim  = 1:8;
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'Speech'
        
        rootdir = fullfile(fn.figs,'ClassSpeech',whichStim,whichClass,['minTrs' num2str(minNtrs)]);
        switch whichClass
            case 'Full'
                rootdir_sum = fullfile(fn.figs,'ClassSpeech',whichStim,'Sum',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'Sum';
            case 'ActVec'
                rootdir_sum = fullfile(fn.figs,'ClassSpeech',whichStim,'SumAV',['minTrs' num2str(minNtrs)]);
                sumclassstr = 'SumAV';
        end
        
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        theseStim  = [4 3 2 1 5 6 7 8];
end

% Load SU data 
q=load(fullfile(rootdir,'CR_each.mat'));
CReach = q.CR;
clear q

% Load subpop results    
% Cat
q=load(fullfile(rootdir,whichCells,['CR_v' whichClass '_' whichCells '.mat']));
CR_Cat = q.CR;
clear q
% Sum
q=load(fullfile(rootdir_sum,whichCells,['CR_v' sumclassstr  '_' whichCells '.mat']));
CR_Sum = q.CR;
clear q


% Also get SU indices and nStim
[~,theseCells] = recallDataParams(whichStim,'each',minNtrs);
nStim          = size(CReach(1,:).Results{:},1);


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen  = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
sqsmall   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)/3*2];


%==========================================================================
%%                         FIG 10 correction
%                    pooled data classifier results
%==========================================================================

% RS indices
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);

% Sort units by peakFR and dprime
[~,iSUdps]  = sort(CReach(iRS,:).dprime,'descend');


%% Collect d' for each stimulus

% - - - RS - - - 
% SU
dp_RS_SU = nan(size(iSUdps,1),nStim);
PC_RS_SU = nan(size(iSUdps,1),nStim);
for ii = 1:numel(iSUdps)
    dp_RS_SU(ii,:) = dp_from_ConfMat(mean(CReach(iRS(iSUdps(ii)),:).Results{:},3,'omitnan'),crAmt);
    PC_RS_SU(ii,:) = diag(mean(CReach(iRS(iSUdps(ii)),:).Results{:},3,'omitnan'));
end

% Cat
dp_RS_Cat  = nan(size(CR_Cat,1),nStim);
PC_RS_Cat  = nan(size(CR_Cat,1),nStim);
sem_Cat    = nan(size(CR_Cat,1),nStim);
NC_Cat     = nan(size(CR_Cat,1),1);
[ncs,incs] = sort(CR_Cat.nC);
for ii = 1:numel(incs)
    dp_RS_Cat(ii,:) = dp_from_ConfMat(mean(CR_Cat(incs(ii),:).Results{:},3,'omitnan'),crAmt);
    PC_RS_Cat(ii,:) = diag(mean(CR_Cat(incs(ii),:).Results{:},3,'omitnan'))';
    sem_Cat(ii,:)   = bootstrap_sem_ConfMat(CR_Cat(incs(ii),:).Results{:},yaxvalue,crAmt);
    NC_Cat(ii)      = ncs(ii);
end

% Sum
dp_RS_Sum  = nan(size(CR_Sum,1),nStim);
PC_RS_Sum  = nan(size(CR_Sum,1),nStim);
sem_Sum    = nan(size(CR_Sum,1),nStim);
NC_Sum     = nan(size(CR_Sum,1),1);
[ncs,incs] = sort(CR_Sum.nC);
for ii = 1:numel(incs)
    dp_RS_Sum(ii,:) = dp_from_ConfMat(mean(CR_Sum(incs(ii),:).Results{:},3,'omitnan'),crAmt);
    PC_RS_Sum(ii,:) = diag(mean(CR_Sum(incs(ii),:).Results{:},3,'omitnan'))';
    sem_Sum(ii,:)   = bootstrap_sem_ConfMat(CR_Sum(incs(ii),:).Results{:},yaxvalue,crAmt);
    NC_Sum(ii)      = ncs(ii);
end

% Correct each to max d' value
dp_RS_SU(dp_RS_SU>dpmaxval)   = dpmaxval;
dp_RS_Cat(dp_RS_Cat>dpmaxval) = dpmaxval;
dp_RS_Sum(dp_RS_Sum>dpmaxval) = dpmaxval;


%% Begin plot: task averages

hf=figure; 
set(hf,'Position',sqsmall)

subplot(2,2,1)
hold on

switch yaxvalue
    case 'dp'
        %~~~~ Add data
        % SU
        dp = CReach(iRS(iSUdps),:).dprime;
        plot(1:length(iSUdps),dp,'.k','MarkerSize',15)
%         plot(1:length(iSUdps),mean(dp_RS_SU,2),'.k','MarkerSize',15)
        
        % Cat
        [ncs,incs] = sort(CR_Cat.nC);
        dp = CR_Cat.dprime(incs);
        dp(dp>dpmaxval) = dpmaxval;
%         plot(ncs,dp,':b','LineWidth',3)
        fill([ncs; flipud(ncs); ncs(1)],[dp - mean(sem_Cat,2); flipud(dp + mean(sem_Cat,2)); dp(1) - mean(sem_Cat(1,:),2)],'b')
        
        % Sum
        [ncs,incs] = sort(CR_Sum.nC);
        dp = CR_Sum.dprime(incs);
        dp(dp>dpmaxval) = dpmaxval;
%         plot(ncs,dp,'-b','LineWidth',3)
        fill([ncs; flipud(ncs); ncs(1)],[dp - mean(sem_Sum,2); flipud(dp + mean(sem_Sum,2)); dp(1) - mean(sem_Sum(1,:),2)],'c')
        
        %~~~~ Finish plot
        ylim([-0.3 dpmaxval])
        ylabel('d'' across stimuli')
        
    case 'PC'
        %~~~~ Add data
        % SU
        plot(1:length(iSUdps),mean(PC_RS_SU,2),'.k','MarkerSize',15)
        % Cat
        fill([ncs; flipud(ncs); ncs(1)],[mean(PC_RS_Cat,2) - mean(sem_Cat,2); flipud(mean(PC_RS_Cat,2) + mean(sem_Cat,2)); mean(PC_RS_Cat(1,:),2) - mean(sem_Cat(1,:),2)],'b')
%         fill([ncs ncs],[mean(PC_RS_Cat,2) - mean(sem_Cat,2) mean(PC_RS_Cat,2) + mean(sem_Cat,2)],':b','LineWidth',3)
        % Sum
%         plot(ncs,mean(PC_RS_Sum,2),'-b','LineWidth',3)
        fill([ncs; flipud(ncs); ncs(1)],[mean(PC_RS_Sum,2) - mean(sem_Sum,2); flipud(mean(PC_RS_Sum,2) + mean(sem_Sum,2)); mean(PC_RS_Sum(1,:),2) - mean(sem_Sum(1,:),2)],'c')
        
        %~~~~ Finish plot
        ylim([0 1])
        ylabel('Average PC')
end

%~~~~ Finish plot
xlabel('Cell N')
grid on
box off
set(gca,'Color','none','xtick',round(linspace(0,length(iSUdps),6)))
xlim([0 numel(iRS)])
title([whichStim ' (SUs ranked by best d'')'])


%% Add results for each stimulus

subplot(2,2,3)
hold on

switch yaxvalue
    case 'dp'
        
        %~~~~ Set stimulus order
        diffs = dp_RS_Sum(end,theseStim) - dp_RS_Cat(end,theseStim);
        [pl_diff,idiff] = sort(diffs,'descend');
        
        %~~~~ Add data
%         plot(1:nStim,dp_RS_Cat(end,theseStim(idiff)),':b','LineWidth',3)
%         plot(1:nStim,dp_RS_Sum(end,theseStim(idiff)),'-b','LineWidth',3)
        plot([1:nStim; 1:nStim],[ dp_RS_Cat(end,theseStim(idiff))-sem_Cat(end,theseStim(idiff)); dp_RS_Cat(end,theseStim(idiff))+sem_Cat(end,theseStim(idiff)) ],'-b','LineWidth',3)
        plot([1:nStim; 1:nStim],[ dp_RS_Sum(end,theseStim(idiff))-sem_Sum(end,theseStim(idiff)); dp_RS_Sum(end,theseStim(idiff))+sem_Sum(end,theseStim(idiff)) ],'.-c','LineWidth',3)
        
        %~~~~  Finish plot
        ylabel('d'' all RS cells')
        ylim([0 dpmaxval])
        
    case 'PC'
        
        %~~~~ Set stimulus order
        diffs = PC_RS_Sum(end,theseStim) - PC_RS_Cat(end,theseStim);
        [pl_diff,idiff] = sort(diffs,'descend');
        
        %~~~~  Add data
%         plot(1:nStim,PC_RS_Cat(end,theseStim(idiff)),':b','LineWidth',3)
%         plot(1:nStim,PC_RS_Sum(end,theseStim(idiff)),'-b','LineWidth',3)
        plot([1:nStim; 1:nStim],[ PC_RS_Cat(end,theseStim(idiff))-sem_Cat(end,theseStim(idiff)); PC_RS_Cat(end,theseStim(idiff))+sem_Cat(end,theseStim(idiff)) ],'-b','LineWidth',3)
        plot([1:nStim; 1:nStim],[ PC_RS_Sum(end,theseStim(idiff))-sem_Sum(end,theseStim(idiff)); PC_RS_Sum(end,theseStim(idiff))+sem_Sum(end,theseStim(idiff)) ],'.-c','LineWidth',3)
        
        %~~~~  Finish plot
        ylabel('PC all RS cells')
        ylim([0 1])
end

%~~~~  Finish plot
grid off
box off
xlabel('Stimulus')
set(gca,'Color','none','xtick',1:nStim,'xticklabel',idiff)
xlim([0.5 nStim+0.5])


% Add plot of Differences between classifiers

subplot(2,2,4)

bar(1:nStim,pl_diff,'b')

switch yaxvalue
    case 'dp'
        ylabel('d'' Sum - d'' Cat')
        ylim([-3.25 1.333])
        % ylim([-1.5 1])
    case 'PC'
        ylabel('PC Sum - PC Cat')
        ylim([-1 0.25])
end

%~~~~  Finish plot
xlabel('Stimulus')
grid off
box off
set(gca,'Color','none','xtick',1:nStim,'xticklabel',idiff)
xlim([0.5 nStim+0.5])

savename = ['F10_' whichStim '_' whichClass '_' num2str(minNtrs) 'trs_' yaxvalue];
print_eps_kp(hf,fullfile(savedir,savename))



end

