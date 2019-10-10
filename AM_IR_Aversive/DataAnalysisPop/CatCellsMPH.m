function CatCellsMPH
% 
% PopMPH
%
% 
%  Intended to help categorize response types.
% 

global AMrates 


%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'z'; %log
Rotate   =    1;
%~~~~~~~~~~~~~~~~~~~~~
clipZ    =    0.5;
%~~~~~~~~~~~~~~~~~~~~~


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen   = [1 scrsz(4) scrsz(3) scrsz(4)];
quartscreen  = [1 scrsz(4) scrsz(3)/2 scrsz(4)/2];

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load PopMPH data
savedir = fullfile(fn.figs,'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

load(fullfile(savedir,'MPHdata.mat'))


%%

switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5 1.75]) ];
end

AllLatencies  = nan(size(ThisData,1),5);
PeakTroughCat = zeros(size(ThisData,1),5);

hf1 = figure;
set(gcf,'Position',fullscreen)

for ir=1:5
    
    data = ThisData(:,1:ceil(1000/AMrates(ir)),ir);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Sort cells
    respdata = [];
    uAbvRnd = ceil(1000/AMrates(ir)/2); %1/4 of period
    
    % do it the slow but easy way
    for iUn=1:size(data,1)
        
        clear uThZ uLat uAbv
        
        if all(isnan(data(iUn,:)))
            respdata = [respdata; nan nan nan ];
        else
            
            %find highest thresh exceeded
            uThZ = Boundaries(find(max(data(iUn,:),[],2) > Boundaries,1,'last'));
            %find time exceeded it
            uLat = find(data(iUn,:) > uThZ,1,'first');
            if mean(data(iUn,1:10))>max(data(iUn,:))  %wrap latency to next period
                uLat = uLat + ceil(1000/AMrates(ir));
            end
%             find ms above
            uAbv = sum(data(iUn,:) > uThZ);
            
            % Save unit data
            respdata = [respdata; uThZ uLat uAbv ];
        end
    end
    [~,i_sorted] = sortrows( respdata, [1 -2 3]);
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Get ready to plot
        
    % ROTATE phase of MPH
    midpoint  = floor(size(data,2)/2);
    if Rotate
        plotdata  = [data(i_sorted,midpoint+1:end) data(i_sorted,1:midpoint) data(i_sorted,midpoint+1:end) data(i_sorted,1:midpoint)]; 
        xtickset  = [0:midpoint:size(plotdata,2)];
        xticklabs = {'-pi' '0' 'pi' '2*pi' 'pi'};
    else 
        plotdata  = [data(i_sorted,:) data(i_sorted,:)];
        xtickset  = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
        xticklabs = [0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))];
    end
    ndp = sum(sum(isnan(plotdata),2)==0);
    
    if clipZ>0
        plotdata(plotdata>clipZ) = clipZ;
        plotdata(plotdata<clipZ) = 0;
    end
    
    % Label NS cells
    flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
    
    % Render subplot
    figure(hf1); hold on
    subplot(1,6,ir+1);
    switch useFR
        case 'z'
            imagesc(plotdata(1:ndp,:))
            caxis([min(Boundaries) max(Boundaries)])
            if clipZ>0
                caxis([0 clipZ])
            end
            cmocean('balance','pivot',0) %curl
        case 'log'
            imagesc(log10(plotdata(1:ndp,:)))
            caxis([0 log10(max(Boundaries))])
            cmocean('gray')
    end
    
    % Add markers to label NS cells
    hold on
    plot(0,find(flagNS),'.g')
    plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.g')
    
    % Finish plot 
    xlim([0 2*ceil(1000/AMrates(ir))])
    ylim([0.5 ndp+0.5])
    set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    set(gca,'xtick',xtickset,'xticklabel',xticklabs)
    
    [r,c]=find(cumsum(respdata(i_sorted(1:ndp),1)==Boundaries)==1);
    set(gca,'ytick',r,'yticklabel',Boundaries)
    
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    
    
    %% Quantify timing of activity
    
    if clipZ>0
        
        % set idx_trough idx_peak
        if Rotate
            idx_tr = floor(midpoint/2 + [0:midpoint-1]);
            idx_pk = floor(idx_tr + midpoint);
        else
            idx_pk = floor(midpoint/2 + [0:midpoint-1]);
            idx_tr = floor(idx_pk + midpoint);
        end
        plot(idx_tr,(ndp+0.5)*ones(size(idx_tr)),'c','LineWidth',5)
        plot(idx_pk,(ndp+0.5)*ones(size(idx_pk)),'m','LineWidth',5)
        
        data_tr  = plotdata(1:ndp,idx_tr);
        data_pk  = plotdata(1:ndp,idx_pk);
        
        count_tr = sum(data_tr>0,2);
        count_pk = sum(data_pk>0,2);
        
        cluD1 = (count_tr-count_pk)./(count_tr+count_pk);
        cluD2 = (count_tr-count_pk);
        
        figure;
        set(gcf,'Position',quartscreen)
        subplot(1,2,1)
        plot(cluD1,cluD2,'ok')
        hold on
        plot([0 0],800/AMrates(ir)*[-1 1],'k--')
        xlabel('(Trough-Peak) / (Trough+Peak)')
        ylabel('Trough - Peak (ms)')
        
        
        % Categorize
        CB = 0;
        TrCat = find(cluD1>CB);
        PkCat = find(cluD1<CB);
        
        PeakTroughCat(i_sorted(TrCat),ir) = -1;
        PeakTroughCat(i_sorted(PkCat),ir) =  1;
        
        
        % Finish subplot 1
        plot(cluD1(PkCat),cluD2(PkCat),'om')
        plot(cluD1(TrCat),cluD2(TrCat),'oc')
        axis square
        set(gca,'Color','none')
        title([num2str(AMrates(ir)) ' Hz'])
        
        % Overlay norm MPHs for each category
        subplot(1,2,2)
        hp1=plot((data(i_sorted(PkCat),:)./max(data(i_sorted(PkCat),:),[],2))','m');
        hold on
        hp2=plot((data(i_sorted(TrCat),:)./max(data(i_sorted(TrCat),:),[],2))','c');
        
        xlim([0 size(data,2)])
        ylim([-1 1].*1.5)
        legend('Peak')
        xlabel(sprintf('Categorized by ms exceeding z=%0.1f\nPref Peak N=%i,  Pref Trough N=%i',clipZ,length(PkCat),length(TrCat)))
        ylabel('Normalized z MPH')
        
        axis square
        set(gca,'Color','none')
        title([num2str(AMrates(ir)) ' Hz'])
        
        
%         savename = sprintf('CatTiming_z%i_%iHz',clipZ*10,AMrates(ir));
%         print_eps_kp(gcf,fullfile(savedir,'CatTiming',savename))
        
        
%         for iu = 1:numel(PkCat)
%             mu = meanphase(find(data(i_sorted(PkCat(iu)),:)>clipZ),size(data,2));
%             AllLatencies(i_sorted(PkCat(iu)),ir) = mu/(2*pi) * size(data,2); %ms
%         end
        for iu = 1:numel(TrCat)
            mu = meanphase(find(data(i_sorted(TrCat(iu)),:)>clipZ),size(data,2));
            AllLatencies(i_sorted(TrCat(iu)),ir) = mu/(2*pi) * size(data,2); %ms
        end
        
    end
    
    
end %ir
% 
% savename = sprintf('PopMPH_%s_Rot%i',useFR,Rotate);
% if clipZ>0
%     savename = [savename '-clipZ' num2str(10*clipZ)];
% end
% 
% print_eps_kp(hf1,fullfile(savedir,savename))


% Add the response category to Units
for iUn = 1:numel(UnitData)
    UnitData(iUn).PeakTroughCat = PeakTroughCat(iUn,:);
end
save(fullfile(fn.processed,'Units_250'),'UnitInfo','UnitData','-v7.3');


% Plot all latencies (of cells that exceed z thresh for that AM rate)
figure; hold on
for iUn = 1:size(AllLatencies,1)
    plot(AllLatencies(iUn,:),1:5,'.-','MarkerSize',20)
end
plot(mean(AllLatencies,1,'omitnan'),1:5,'-k','LineWidth',8)

set(gca,'ytick',1:5,'yticklabel',AMrates,'Color','none')
ylabel('AM rate')
xlabel('Mean phase (ms)')
title('Trough Cells')

% print_eps_kp(gcf,fullfile(savedir,['LatShift_MeanPhase_clipZ_Peak' num2str(clipZ*10)]))


% Load RCorr results
% q=load(fullfile(fn.figs,'RCorr','exclOnset','PCMat_10trTemp_b.mat'));
% PCMat = q.PCMat;
% clear q
% if size(PCMat,3) ~= numel(UnitData)
%     warning('cant sort by RCorr accuracy')
%     keyboard
% end




end