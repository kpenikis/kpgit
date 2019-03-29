function plotMPH_aggregatePopulations(SUBJECT,SESSION)
% 
% Must be single session.
% Flips through trials for specified stimulus. Plots spike times of
% selected RS units and all NS units.
% 
% KP, 2019-03
% 


close all
fn = set_paths_directories(SUBJECT,SESSION,1);


% Get raster data
[StimResp, UnitData, UnitInfo, Info, TrialData] = collectRasterDataSession(SUBJECT,SESSION);

% SELECT UNITS
% Sequence of response peaks, by eye, from Session raster
% %
cluIDs = [1650 1139 771 1230]; %68
% %
[~,~,UIidxRS] = intersect(cluIDs,UnitInfo.Clu,'stable');
% Also get indices of NS units
UIidxNS = find(strcmp(UnitInfo.SpkShape,'NS'));






%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
tallrect  = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];

UnitColors = cmocean('phase',length(cluIDs)+3);

        
% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


%% Select stimuli to plot

for stid = 3
    
    % Prepare the figure
    hf(stid) = figure;
    set(hf(stid),'Position',tallrect.* [1 1 figwidthscales(stid) 1],'NextPlot','add')
    hold on
%     hs(1)=subplot(9,1,1);
%     set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[],'Color','none');
%     box off
    hs(2)=subplot(9,1,1:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'ytick',[],...
        'xticklabel',[0 (figwidthscales(stid)*1000)],'Color','none','tickdir','out');
    xlabel('Time (ms)')
    box on
    
    
    % Plot rasters/PSTHs
    
    for it =  1:size(StimResp(stid).raster,1)
        
        N_un = 0;
        
        pd_x = 0:(1000/AMrates(stid-1)):figwidthscales(stid)*1000;
        pd_x = [pd_x; pd_x];
        pd_y = [zeros(1,size(pd_x,2))+0.5; (numel(UIidxRS)+numel(UIidxNS)+1).*ones(1,size(pd_x,2))];
        
        subplot(hs(2)); cla;
        plot(pd_x,pd_y,'k-')
        hold on
        
        % RS units
        for iClu = UIidxRS'
            
            N_un = N_un + 1;
            
            % Convert spiketime format
            raster_x=[]; raster_y=[];
            raster_x = find(StimResp(stid).raster(it,:,iClu));
            raster_y = N_un * ones(1,sum(StimResp(stid).raster(it,:,iClu)));
            
            % Add to plot
            plot([raster_x; raster_x], [raster_y; raster_y] + 0.45.*[-ones(size(raster_y)); ones(size(raster_y))],...
                'Color',UnitColors(N_un+1,:),'LineWidth',5)
            
        end %iClu
        
        % NS units
        for iClu = UIidxNS'
            
            N_un = N_un + 1;
            
            patch([0 figwidthscales(stid)*1000 figwidthscales(stid)*1000 0],...
                [N_un-0.5 N_un-0.5 N_un+0.5 N_un+0.5],...
                0.2*[1 1 1],'EdgeColor','none','FaceAlpha',0.1)
            
            % Convert spiketime format
            raster_x=[]; raster_y=[];
            raster_x = find(StimResp(stid).raster(it,:,iClu));
            raster_y = N_un * ones(1,sum(StimResp(stid).raster(it,:,iClu)));
            
            % Add to plot
            plot([raster_x; raster_x], [raster_y; raster_y] + 0.45.*[-ones(size(raster_y)); ones(size(raster_y))],...
                'Color',UnitColors(N_un+1,:),'LineWidth',5)
            
        end %iClu
        
        set(gca,'ylim',[0.5 (numel(UIidxRS)+numel(UIidxNS)+0.5)])
        pause(1.5)
        
    end %it
    
end %stid



end