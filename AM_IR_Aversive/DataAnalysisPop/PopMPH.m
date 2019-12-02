function PopMPH
% 
% PopMPH
%
%  View MPHs of all cells in population.
% 

global AMrates useFR Boundaries


%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
Rotate   =    1;
%~~~~~~~~~~~~~~~~~~~~~
clipZ    =    0;
%~~~~~~~~~~~~~~~~~~~~~
plMPH    =    1;
spks     =    0;
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
set(0,'DefaultAxesFontSize',8)

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


%% Filter data to just one session

% keyboard
% SUBJ = 'AAB_265054';
% SESS = 'Apr02-AM';
% 
% theseUnits = find(strcmp(UnitInfo.Subject,SUBJ) & strcmp(UnitInfo.Session,SESS));
% 
% ThisData = ThisData(theseUnits,:,:);
% FR_Warn  = FR_Warn(theseUnits,:);


%%

if plMPH
    
    hf1 = figure;
    set(gcf,'Position',fullscreen)
    hold on
    
%     keyboard %check how you want to sort
%     [i_sorted,sortdata] = sort_thLat(ThisData(:,1:250,2));
    
    % Warn
    subplot(1,6,1);
    plotWarnPopH
    
    % Each AM rate
    for ir=1:5
        
%         data = ThisData(:,1:ceil(1000/AMrates(ir)),ir);
%         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         % Sort cells
%         [i_sorted,sortdata] = sort_thLat(data);

        subplot(1,6,ir+1);
        plotThisMPH  %_2pd
    end
    
    % Save figure
    savename = sprintf('PopMPH_%s_noshift',useFR);
    % savename = sprintf('PopMPH_%s_Rot%i',useFR,Rotate);
    if clipZ>0
        savename = [savename '-clipZ' num2str(10*clipZ)];
    end
    
    print_eps_kp(gcf,fullfile(savedir,savename))
    set(gcf,'PaperOrientation','landscape')
    print(gcf,'-dpdf','-r500','-fillpage', fullfile(savedir,savename))
    
end %if plMPH


return

%%
keyboard
if spks

hf4 = figure;
set(gcf,'Position',fullscreen)

for itr = 1:size(sp_trs,4)

    clf
    hold on
    
    for ir=1:5
        
        data   = ThisData(:,1:ceil(1000/AMrates(ir)),ir);
        spdata = sp_trs(:,1:ceil(1000/AMrates(ir)),ir,itr);
        
        % Sort cells by Latency
        [i_sorted,sortdata] = sort_thLat(data);
        
        
        % ROTATE phase of MPH
        if Rotate
            midpoint  = floor(size(data,2)/2);
            plotdata  = [data(i_sorted,midpoint+1:end)   data(i_sorted,1:midpoint)   data(i_sorted,midpoint+1:end)   data(i_sorted,1:midpoint)];
            plotspks  = [spdata(i_sorted,midpoint+1:end) spdata(i_sorted,1:midpoint) spdata(i_sorted,midpoint+1:end) spdata(i_sorted,1:midpoint)];
            xtickset  = [0:midpoint:size(plotdata,2)];
            xticklabs = {'-pi' '0' 'pi' '2*pi' 'pi'};
        else
            plotdata  = [data(i_sorted,:)   data(i_sorted,:)];
            plotspks  = [spdata(i_sorted,:) spdata(i_sorted,:)];
            xtickset  = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
            xticklabs = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
        end
        
        ndp = sum(sum(isnan(plotdata),2)==0);
        
        
        % Label NS
        flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
        
        % Get corresponding spike data
        
        % Now plot raster
        subplot(1,6,ir+1);
        
        [x,y] = find(plotspks(1:ndp,:)');
        plot([x x]',y'+[-0.5; 0.5],'-k')
        set(gca,'ydir','reverse')
        
        % Add markers to label NS cells
        hold on
        plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
        plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.','Color',[0.01 0.57 0.44])
        
        % Finish plot
        xlim([0 2*ceil(1000/AMrates(ir))])
        ylim([0.5 ndp+0.5])
        set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
        set(gca,'xtick',xtickset,'xticklabel',xticklabs)
        
        [r,~]=find(cumsum(sortdata(i_sorted(1:ndp),1)==Boundaries)==1);
        set(gca,'ytick',r,'yticklabel',Boundaries)
        
        title([num2str(AMrates(ir)) ' Hz'])
        box off
        axis fill
        
        
        % Plot average histogram
        subplot(1,6,1);
        hold on 
        plot(cumsum(sum(plotspks(1:ndp,1:size(plotspks,2)/2),1),2),'Color',colors(ir+1,:),'LineWidth',2)
%         axis square
        
%         plot(1:size(frdata,2),ir+0.01.*mean(spdata(1:ndp,:),1).*1000,'k','LineWidth',2)
%         plot(1:size(frdata,2),ir+0.2.*median(frdata,1,'omitnan'),'b','LineWidth',2)
        
    end %ir
    
    xlim([0 500])
    ylim([0 800])
    ylabel('Cumulative spike count')
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    set(gca,'xtick',[31 63 125 250 500])
    xlabel('Time (ms)')
    title(['random trial ' num2str(itr)])
    
    set(gcf,'PaperOrientation','landscape')
%     print_eps_kp(gcf,fullfile(savedir,['PopMPH_log_Spks_Tr-' num2str(itr) ]))
    print(gcf,'-dpdf','-r500','-fillpage', fullfile(savedir,'Trials',sprintf('PopMPH_%s_Spks_Tr-%i',useFR,itr)) )
    
end %itr

end %if spks



%% plot 6 trials of one rate

keyboard

convwin     = 10;
convwin     = ones(1,convwin).*(1/convwin);

hf4 = figure;
set(gcf,'Position',fullscreen)

ir=2;

for itr = 1:6
    
    data   = ThisData(:,1:ceil(1000/AMrates(ir)),ir);
    spdata = sp_trs(:,1:ceil(1000/AMrates(ir)),ir,itr);
    
    % Sort cells by Latency
    [i_sorted,sortdata] = sort_thLat(data);
    
    
    % ROTATE phase of MPH
    if Rotate
        midpoint  = floor(size(data,2)/2);
        plotdata  = [data(i_sorted,midpoint+1:end)   data(i_sorted,1:midpoint)   data(i_sorted,midpoint+1:end)   data(i_sorted,1:midpoint)];
        plotspks  = [spdata(i_sorted,midpoint+1:end) spdata(i_sorted,1:midpoint) spdata(i_sorted,midpoint+1:end) spdata(i_sorted,1:midpoint)];
        xtickset  = [0:midpoint:size(plotdata,2)];
        xticklabs = {'-pi' '0' 'pi' '2*pi' 'pi'};
    else
        plotdata  = [data(i_sorted,:)   data(i_sorted,:)];
        plotspks  = [spdata(i_sorted,:) spdata(i_sorted,:)];
        xtickset  = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
        xticklabs = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
    end
    
    ndp = sum(sum(isnan(plotdata),2)==0);
    
    
    % Get corresponding spike data
    
    % Now plot raster
    subplot(4,6,itr);
    
    % psth
    plot(sum(plotspks(1:ndp,:),1),'k')
    smFR = smoothFR(sum(plotspks(1:ndp,:),1)./ndp,20);
    hold on
    plot(smFR,'b')
    foo = conv(sum(plotspks(1:ndp,:),1),convwin);
    ps_conv = foo(floor(numel(convwin)/2)+(0:size(plotspks,2)-1));
    plot(midpoint+(0:ceil(1000/AMrates(ir))),ps_conv(midpoint+(0:ceil(1000/AMrates(ir)))),'m','LineWidth',2)
    
    % raster
%     [x,y] = find(plotspks(1:ndp,:)');
%     plot([x x]',y'+[-0.5; 0.5],'-k')
%     set(gca,'ydir','reverse')
    
    % Add markers to label NS cells
    %         flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
    %         hold on
    %         plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
    %         plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.','Color',[0.01 0.57 0.44])
    
    % Finish plot
    xlim([0 2*ceil(1000/AMrates(ir))])
%     ylim([0.5 ndp+0.5])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    set(gca,'xtick',xtickset,'xticklabel',xticklabs)
    ylim([0 10])
%     [r,~]=find(cumsum(sortdata(i_sorted(1:ndp),1)==Boundaries)==1);
%     set(gca,'ytick',r,'yticklabel',Boundaries)
    
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    
    set(gcf,'PaperOrientation','landscape')
    
    
end %itr

savename = sprintf('PopMPH_%sPSTH_%ihz_Tr-1-%i',useFR,AMrates(ir),itr);
print_eps_kp(gcf,fullfile(savedir,'Trials',savename))
print(gcf,'-dpdf','-r500','-fillpage', fullfile(savedir,'Trials',savename) )


end