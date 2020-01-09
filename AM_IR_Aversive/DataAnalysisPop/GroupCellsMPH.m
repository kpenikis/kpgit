function GroupCellsMPH
% 
% GroupCellsMPH
%
%  Group cells in population by their MPHs. 
% 

global AMrates useFR Boundaries


%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
useGrp   =   'peakFR'; 'dynRange'; 'tuning'; 'phase'; 
%~~~~~~~~~~~~~~~~~~~~~
ngrps    = 5;

close all

%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load PopMPH data
savedir = fullfile([fn.figs '_v1'],'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

load(fullfile(savedir,'MPHdata.mat'))

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end


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


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];


%%

% Reshape data / concatenate MPHs 

foo = reshape(ThisData,[size(ThisData,1) size(ThisData,2)*size(ThisData,3) 1]);

CatData = nan(size(ThisData,1),sum(ceil(1000./AMrates)));

for iUn = 1:size(ThisData,1)
    
    foo_abbr = foo(iUn,~isnan(foo(iUn,:)));
    
    if length(foo_abbr) == sum(ceil(1000./AMrates))
        CatData(iUn,:) = foo_abbr;
    end
end


% Pull numbers from each MPH
DynRange = nan(size(ThisData,1),size(ThisData,3),2);
MeanRate = nan(size(ThisData,1),size(ThisData,3));
PeakLat  = nan(size(ThisData,1),size(ThisData,3));

for iUn = 1:size(ThisData,1)
    for ist = 1:size(ThisData,3)
        
        DynRange(iUn,ist,1) = min(ThisData(iUn,:,ist));
        [DynRange(iUn,ist,2),PeakLat(iUn,ist)] = max(ThisData(iUn,:,ist));
        MeanRate(iUn,ist)   = mean(ThisData(iUn,:,ist),'omitnan');
        
    end
end


%%

switch useGrp
    
    %~~~~~~~~~~~~~~~~~~~
    case 'peakFR'
    % Group by peak rate
        
        quant_bounds_p = linspace(0,1,ngrps+1);
        
        GrpData = max(DynRange(:,:,2),[],2);
        
        GrpBounds = quantile(GrpData,quant_bounds_p);
        
        figure;
        plot(sort(GrpData),'.b')
        hold on
        plot([0 size(ThisData,1)],[GrpBounds' GrpBounds'],'k')
        xlabel('Unit #')
        ylabel('Peak Rate')
        title('Rough grouping by peak FR')
        
        
    %~~~~~~~~~~~~~~~~~~~
    case 'dynRange'
    % Group by dynamic range
        
        quant_bounds_p = linspace(0,1,ngrps+1);
        
        GrpData = mean( diff(DynRange,1,3) ,2);
        
        GrpBounds = quantile(GrpData,quant_bounds_p);
        
        figure;
        plot(sort(GrpData),'.c')
        hold on
        plot([0 size(ThisData,1)],[GrpBounds' GrpBounds'],'k')
        xlabel('Unit #')
        ylabel('Mean dynamic range')
        title('Rough grouping by FR dynamic range')
        
        
    %~~~~~~~~~~~~~~~~~~~
    case 'tuning'
    % Group by dynamic range of tuning curve
        
        quant_bounds_p = linspace(0,1,ngrps+1);
        
        GrpData = (max(MeanRate,[],2) - min(MeanRate,[],2)) ./ ([UnitData.BaseFR]'+0.001);
        
        GrpBounds = quantile(GrpData,quant_bounds_p);
        
        figure;
        plot(sort(GrpData),'.b')
        hold on
        plot([0 size(ThisData,1)],[GrpBounds' GrpBounds'],'k')
        xlabel('Unit #')
        ylabel('Peak Rate')
        title('Rough grouping by tuning range')
        
        
	%~~~~~~~~~~~~~~~~~~~
    case 'phase'
    % Group by avg phase 2 to 8 Hz
        
        quant_bounds_p = linspace(0,1,ngrps+1);
        
        allPhases = nan(numel(UnitData),5);
        for iUn = 1:numel(UnitData)
            allPhases(iUn,:) = UnitData(iUn).Phase_spk(2:6);
        end
        GrpData = circ_mean(deg2rad(allPhases(:,1:3))-pi,[],2)+pi;
        
        GrpBounds = quantile(GrpData,quant_bounds_p);
        
        figure;
        plot(sort(GrpData),'.r')
        hold on
        plot([0 size(ThisData,1)],[GrpBounds' GrpBounds'],'k')
        xlabel('Unit #')
        ylabel('Mean phase (2,4,8 Hz)')
        title('Rough grouping by preferred phase')
        
end


% Label cells with group #
GroupLabel = nan(size(ThisData,1),1);
for iq = 2:numel(GrpBounds)
    GroupLabel( GrpData>GrpBounds(iq-1) & GrpData<=GrpBounds(iq) ) = iq-1;
end


% Plot MPHs of each group

hf = figure;
set(gcf,'Position',fullscreen)
hold on
    
for iGrp = 1:ngrps
    
    iUns = find(GroupLabel==iGrp);
    
    data = CatData(iUns,:);
    [i_sorted,sortdata] = sort_thLat(data);
    
    subplot(1,ngrps,iGrp);
    hold on
    
    plotdata  = data(i_sorted,:) ;
    xtickset  = [0 cumsum(ceil(1000./AMrates))];
    xticklabs = [0 cumsum(ceil(1000./AMrates))];
    
    ndp = sum(sum(isnan(plotdata),2)==0);
    
    % Render plot
    switch useFR
        case 'z'
            imagesc(plotdata(1:ndp,:))
            caxis([min(Boundaries) max(Boundaries)])
            cmocean('balance','pivot',0) %curl
        case 'log'
            imagesc(log10(plotdata(1:ndp,:)))
            %         caxis([0 log10(max(Boundaries))+0.25])
            caxis([0 1.75])
            cmocean('gray')
    end
    
    % Add markers to label NS cells
    flagNS = UnitInfo.TroughPeak(iUns(i_sorted(1:ndp)))<0.5;
    hold on
    plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
    plot(size(plotdata,2),find(flagNS),'.','Color',[0.01 0.57 0.44])
    
    % Finish plot
    xlim([0 size(plotdata,2)])
    ylim([0.5 ndp+0.5])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    set(gca,'xtick',xtickset,'xticklabel',xticklabs)
    
    BoundMat = cumsum(sortdata(i_sorted(1:ndp),1)>=Boundaries);
    ytplc = []; ytlab = [];
    for ib = numel(Boundaries):-1:1
        yUn = find(BoundMat(:,ib)==1,1,'first');
        if ~ismember(yUn,ytplc)
            ytplc = [ytplc yUn];
            ytlab = [ytlab Boundaries(ib)];
        end
    end
    set(gca,'ytick',fliplr(ytplc),'yticklabel',fliplr(ytlab),'ydir','reverse')
    
    % CluLabels = [UnitData(theseUnits(i_sorted)).Clu];
    % set(gca,'ytick',theseUnits,'yticklabel',CluLabels)
    
    title([' Group ' num2str(iGrp) ' MPHs'])
    box off
    axis fill
    
end


% Save figure
savedir = fullfile(fn.figs,'PopMPH','Grouped');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
savename = sprintf('PopMPH_%s_%s_%i',useFR,useGrp,ngrps);

print_eps_kp(gcf,fullfile(savedir,savename))
% set(gcf,'PaperOrientation','landscape')
% print(gcf,'-dpdf','-r500','-fillpage', fullfile(savedir,savename))


% Save Group Labels
savename = sprintf('GroupLabels_%s_%i',useGrp,ngrps);
save(fullfile(savedir,savename),'GroupLabel','-v7.3');




end

