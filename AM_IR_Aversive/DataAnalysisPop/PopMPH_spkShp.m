function PopMPH_spkShp
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
whichCells = 'NS'; 'AAB_265054'; 'WWWlf_253395'; 'AAB_265058'; 'WWWf_253400'; 'Jan25'; 'all'; 
%~~~~~~~~~~~~~~~~~~~~~
sortHere = 1;
sortBy   = 'tMax4'; 'baseFR'; 'FF4'; 
%~~~~~~~~~~~~~~~~~~~~~
labelNS  = 0;
%~~~~~~~~~~~~~~~~~~~~~

% close all

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
q=load(fullfile(fullfile(fn.figs,'StimClass'),'Cell_Time_Trial_Stim_simtrs'));
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;


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


% Label NS cells
flagNS = UnitInfo.TroughPeak<0.43;

switch whichCells
    case 'NS'
        theseCells = find(flagNS);
    case 'RS'
        theseCells = find(~flagNS);
    case 'all'
        theseCells = 1:size(ThisData,1);
    case {'Apr02' 'Apr07' 'Apr09' 'Apr11' 'Apr15' 'Mar28' 'Mar30' 'Jan17' 'Jan21' 'Jan25'} 
        theseCells = find(strcmp({UnitData.Session},[whichCells '-AM']));
    case {'AAB_265054' 'WWWf_253400' 'AAB_265058' 'WWWlf_253395'}
        theseCells = find(strcmp({UnitData.Subject},whichCells));
end

convwin = 10;
convwin = ones(1,convwin).*(1/convwin);


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


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
    savedir = fullfile(fn.figs,'PopMPH','Dec2019');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    savename = sprintf('PopMPH_%s_%s_%s',useFR,whichCells,sortBy);
    if clipZ>0
        savename = [savename '-clipZ' num2str(10*clipZ)];
    end
    
%     print_eps_kp(gcf,fullfile(savedir,savename))
    set(hf1,'PaperOrientation','landscape')
    print(hf1,'-dpdf','-r500','-fillpage', fullfile(savedir,savename))
    
end %if plMPH


return

%%

if spks

hf4 = figure;
set(gcf,'Position',fullscreen)

for itr = 1:size(sp_trs,4)

    clf
    hold on
    
    % Each AM rate
    for ir=1:5
        
        %
        data = ThisData(:,1:ceil(1000/AMrates(ir)),ir);
        
        % Length of period
        spdata = sp_trs(:,1:ceil(1000/AMrates(ir)),ir,itr);
        
        % Full 500 ms duration (different trials than raster plots)        
%         spdata = Cell_Time_Trial_Stim(:,:,itr,ir+1);
        
        
        % %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % % Sort cells
        [i_sorted,sortdata] = sort_thLat(data);
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Get ready to plot
        
        % ROTATE phase of MPH
        % if Rotate
        %     midpoint  = floor(size(data,2)/2);
        %     plotdata  = [data(i_sorted,midpoint+1:end) data(i_sorted,1:midpoint) data(i_sorted,midpoint+1:end) data(i_sorted,1:midpoint)];
        %     xtickset  = [0:midpoint:size(plotdata,2)];
        %     xticklabs = {'-pi' '0' 'pi' '2*pi' 'pi'};
        % else
        plotdata  = data(i_sorted,:);
        plotspks  = spdata(i_sorted,:);
        xtickset  = [0 ceil(1000/AMrates(ir))/2 ceil(1000/AMrates(ir))];
        xticklabs = {'0' 'pi' '2*pi'};
        % end
        ndp = sum(sum(isnan(plotdata),2)==0);
        
        
        % Label NS cells
        flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
        
        
        % Now plot raster
        subplot(1,6,ir+1);
        
        [x,y] = find(plotspks(1:ndp,:)');
        plot([x x]',y'+[-0.5; 0.5],'-k','LineWidth',2)
        set(gca,'ydir','reverse')
        
        
%         % Render plot
%         switch useFR
%             case 'z'
%                 imagesc(plotdata(1:ndp,:))
%                 caxis([min(Boundaries) max(Boundaries)])
%                 cmocean('balance','pivot',0) %curl
%             case 'log'
%                 imagesc(log10(plotdata(1:ndp,:)))
%                 %         caxis([0 log10(max(Boundaries))+0.25])
%                 caxis([0 1.75])
%                 cmocean('-gray')
%         end
        
        % Add markers to label NS cells
        hold on
        plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
        plot(ceil(1000/AMrates(ir)),find(flagNS),'.','Color',[0.01 0.57 0.44])
        
        % Finish plot
        xlim([0 ceil(1000/AMrates(ir))])
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
        set(gca,'ytick',fliplr(ytplc),'yticklabel',fliplr(ytlab))
        
        % CluLabels = [UnitData(theseUnits(i_sorted)).Clu];
        % set(gca,'ytick',theseUnits,'yticklabel',CluLabels)
        
        title([num2str(AMrates(ir)) ' Hz'])
        box off
        axis fill
        
        
        % Plot average histogram
        subplot(1,6,1);
        plot([0 size(plotspks,2)],[ir ir],':k')
        hold on 
        convplot = conv([mean(plotspks(end-length(convwin)+1:end),1,'omitnan') mean(plotspks,1,'omitnan') mean(plotspks(1:length(convwin)),1,'omitnan')],convwin,'same');
        plot(convplot(length(convwin)+(1:size(plotspks,2)))*30+ir,'k','LineWidth',2)
        
        
    end %ir
    
    xlim([0 500])
    ylim([0.5 5.5])
%     ylabel('Cumulative spike count')
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    set(gca,'xtick',[0 500],'ytick',[1:5],'yticklabel',AMrates)
    xlabel('Time (ms)')
    title(['random trial ' num2str(itr)])
    box off
    
    set(gcf,'PaperOrientation','landscape')
    pause(0.5)
%     print_eps_kp(gcf,fullfile(savedir,['PopMPH_log_Spks_Tr-' num2str(itr) ]))
%     print(gcf,'-dpdf','-r500','-fillpage', fullfile(savedir,'Trials',sprintf('PopMPH_%s_Spks_Tr-%i',useFR,itr)) )
    
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
    
    plotdata  = [data(i_sorted,:)   data(i_sorted,:)];
    plotspks  = [spdata(i_sorted,:) spdata(i_sorted,:)];
    xtickset  = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
    xticklabs = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
    
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