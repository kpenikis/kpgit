function PopMPH_bySession
% 
% PopMPH
%
%  View MPHs of all cells in population.
% 

global AMrates useFR Boundaries


%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
%~~~~~~~~~~~~~~~~~~~~~
clipZ    =    0;
%~~~~~~~~~~~~~~~~~~~~~
plMPH    =    1;
spks     =    0;
%~~~~~~~~~~~~~~~~~~~~~
sortHere = 1;
sortBy   = 'tMax4'; 'baseFR'; 'FF4'; 
%~~~~~~~~~~~~~~~~~~~~~
labelNS  = 0;
%~~~~~~~~~~~~~~~~~~~~~

close all

%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load PopMPH data 
%  (created in savePopMPHdata)
savedir = fullfile(fn.figs,'PopMPH'); 
if ~exist(savedir,'dir')
    mkdir(savedir)
end

load(fullfile(savedir,'MPHdata_noshift.mat'))

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
% q=load(fullfile(fullfile(fn.figs,'StimClass'),'Cell_Time_Trial_Stim_simtrs'));
% Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
% if size(Cell_Time_Trial_Stim,1)==257 && ~exist('Un_Indices','var')
%     Un_Indices = 1:257;
% end


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


%% Data settings

switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
%         Boundaries  = [-1 round(10.^[0.5]) ];
end

convwin = 10;
convwin = ones(1,convwin).*(1/convwin);


%% Sort sessions by number of cells
[Sessions,~,idxSess] = unique(UnitInfo(:,1:2),'stable');

y = histcounts(idxSess);

[NunS,iss] = sort(y,'ascend');

% Only sessions with >=10 units
startSess = find(NunS>10,1,'first');
nCells = sum(NunS(startSess:end));


%% Preallocate
PopPlotWarn=[];
PopPlotAM = { nan(nCells,500) nan(nCells,250) nan(nCells,125) nan(nCells,62) nan(nCells,31) };

for ii = startSess:numel(iss)
    
    theseCells = find(idxSess==iss(ii));
    
    clear i_sorted
    if sortHere
        %         [i_sorted,sortdata] = sort_thLat(ThisData(theseCells,:,2));
        switch sortBy
            case 'baseFR'
                [sortdata,i_sorted] = sort([UnitData(theseCells).BaseFR]);
            case 'tMax4'
                [pk,ipk] = max(ThisData(theseCells,1:250,2),[],2);
                [sortdata,i_sorted] = sortrows([ipk pk],[-1 2]);
            case 'FF4'
                FF = var(permute(sum(Cell_Time_Trial_Stim(theseCells,1:250,:,3),2),[1 3 2 4]),[],2,'omitnan') ...
                    ./ mean(permute(sum(Cell_Time_Trial_Stim(theseCells,1:250,:,3),2),[1 3 2 4]),2,'omitnan');
                [sortdata,i_sorted] = sort(FF);
        end
    end
    
    
    %% Each AM rate
    
    for ir=1:5
        
        plotdata = ThisData(theseCells(i_sorted),1:floor(1000/AMrates(ir)),ir);
        
        if clipZ>0
            plotdata(plotdata>clipZ) = clipZ;
            plotdata(plotdata<clipZ) = 0;
        end
        
        % Add to population plot data
        PopPlotAM{ir}(size(PopPlotWarn,1)+(1:size(plotdata,1)),:) = plotdata;
        
    end
    
    
    %% Warn
    switch useFR
        case 'z'
            plotdata = zFR_Warn(theseCells(i_sorted),:);
        case 'log'
            plotdata = FR_Warn(theseCells(i_sorted),:);
    end
    
    if clipZ>0
        plotdata(plotdata>clipZ) = clipZ;
        plotdata(plotdata<clipZ) = 0;
    end
    
    % Add to population plot data
    PopPlotWarn = [PopPlotWarn; plotdata];
    
    
end %step through sessions


%% Plot now

hf1 = figure;
set(gcf,'Position',fullscreen)
hold on

% Warn
subplot(1,6,1);
switch useFR
    case 'z'
        imagesc(PopPlotWarn)%(~any(isnan(PopPlotWarn),2),:))
        caxis([min(Boundaries) max(Boundaries)])
        cmocean('balance','pivot',0) %curl
    case 'log'
        imagesc(log10(PopPlotWarn))
%         caxis([0 log10(max(Boundaries))+0.25])
        caxis([0 1.75])
        cmocean('-gray')
end

% Finish plot
xlim([0 size(PopPlotWarn,2)])
ylim([0.5 size(PopPlotWarn,1)+0.5])
set(gca,'tickdir','out','ticklength',[0.01 0.01],'Color','none')
set(gca,'xtick',[0 size(PopPlotWarn,2)])
set(gca,'ytick',cumsum(NunS(startSess:end)),'yticklabel', Sessions.Session(iss(startSess:end))')
ylabel(sortBy)


% AM rates
for ir=1:5
    
    subplot(1,6,ir+1);
    
    switch useFR
        case 'z'
            imagesc(PopPlotAM{ir})%(~any(isnan(PopPlotWarn),2),:))
            caxis([min(Boundaries) max(Boundaries)])
            cmocean('balance','pivot',0) %curl
        case 'log'
            imagesc(log10(PopPlotAM{ir}))
            %         caxis([0 log10(max(Boundaries))+0.25])
            caxis([0 1.75])
            cmocean('-gray')
    end
    
    % Finish plot
    xlim([0 size(PopPlotAM{ir},2)])
    ylim([0.5 size(PopPlotAM{ir},1)+0.5])
    set(gca,'tickdir','out','ticklength',[0.01 0.01],'Color','none')
    set(gca,'xtick',[0 size(PopPlotAM{ir},2)])
    set(gca,'ytick',cumsum(NunS(startSess:end)),'yticklabel',[])
    
end


%% Save figure
savedir = fullfile(fn.figs,'PopMPH','bySess');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('PopMPH_bySess_%s_%s_abv10',useFR,sortBy);
if clipZ>0
    savename = [savename '-clipZ' num2str(10*clipZ)];
end

%     print_eps_kp(gcf,fullfile(savedir,savename))
    set(hf1,'PaperOrientation','landscape')
    print(hf1,'-dpdf','-r500','-fillpage', fullfile(savedir,savename))


end