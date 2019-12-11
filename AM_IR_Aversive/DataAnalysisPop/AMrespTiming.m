function AMrespTiming
% 
% AMrespTiming
%
%  Group cells in population by their MPHs. 
% 

global AMrates useFR Boundaries

%~~~~~~~~~~~~~~~~~~~~~
plotF1   = 0;
plotF2   = 0;
plotF3   = 1;
plotF4   = 0;
plotF5   = 0;
%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
%~~~~~~~~~~~~~~~~~~~~~
whichLat = 'Mid'; 'Offset'; 'Onset'; 
%~~~~~~~~~~~~~~~~~~~~~

AMrates = [2 4 8 16 32];

% close all


%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

%% CellTypes

CellTypes = {'RS' 'NS'};
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<0.43 & [UnitData.BaseFR]'>3);


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
sz1  = [1 scrsz(4)/4*3 scrsz(3)/2 scrsz(4)/4*3];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [ 181   0  52;...
           255  87  51;...
           255 153   0]./255;
       
savedir = fullfile(fn.figs,'RespChar','CellType');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
% 3rd dim: [mu VS RS]
PhaseData = nan(numel(UnitData),5,3);

for iUn = 1:numel(UnitData)
    PhaseData(iUn,:,1) = UnitData(iUn).Phase_spk(2:6);
    PhaseData(iUn,:,2) = UnitData(iUn).VSdata_spk(1,2:6);
    PhaseData(iUn,:,3) = UnitData(iUn).VSdata_spk(2,2:6);
end


%% N cells synchronized

if plotF1
        
    % Sort data
    sortdata = PhaseData(:,:,3)>13.1;
    [sortdata,i_sorted] = sortrows( sortdata, -1*(1:5) );
    
    % Plot
    figure;
    set(gcf,'Position',sz1)
    
    subplot(1,5,1:3)
    imagesc(PhaseData(i_sorted,:,3)>13.1)
    cmocean('-grey')
    set(gca,'ydir','normal','ylim',[0 numel(i_sorted)])
    set(gca,'xtick',1:5,'xticklabel',AMrates)
    xlabel('AM rate')
    ylabel('Count cells')
    
    
    % Plot FRs next to cells, to get a sense of how much each contributes
    
    AMFR = nan(1,numel(i_sorted));
    for iUn=1:numel(i_sorted)
        AMFR(iUn) = mean(mean(UnitData(iUn).FR_raw_tr(:,7:end),1,'omitnan'),'omitnan');
    end
    
    BaseSubFR = AMFR(i_sorted)-[UnitData(i_sorted).BaseFR];
    
    
    subplot(1,5,4)
    cm=cmocean('speed',60);
    scatter(AMFR(i_sorted),1:numel(i_sorted),20,cm(round(AMFR(i_sorted))+1,:),'filled')
    ylim([0 numel(i_sorted)])
    xlim([0 60])
    set(gca,'ytick',[],'Color','none')
    xlabel('Irr FR')
    
    subplot(1,5,5)
    cm=cmocean('-balance',50);
    scatter(BaseSubFR,1:numel(i_sorted),20,cm(round(BaseSubFR)+25,:),'filled')
    ylim([0 numel(i_sorted)])
    xlim([-25 25])
    set(gca,'ytick',[],'Color','none')
    xlabel('Base sub FR')
    
    
    % Save
    savename = 'PopAMsync_wFR';
    % print_eps_kp(gcf,fullfile(savedir,savename))
    set(gcf,'PaperOrientation','portrait')
    print(gcf,'-dpdf','-bestfit','-r500', fullfile(savedir,savename) )
    
end


%% Phase by VS/nSpks
if plotF2
    
    radvector = 0:15:360;
    
    hfsum=figure; hold on
    set(gcf,'Position',fullscreen)
    
    for ir = 1:size(PhaseData,2)
        
        Phs = PhaseData(PhaseData(:,ir,3)>13.1,ir,1);
        VSs = PhaseData(PhaseData(:,ir,3)>13.1,ir,2);
        RSs = PhaseData(PhaseData(:,ir,3)>13.1,ir,3);
%         Phs = PhaseData(:,ir,1);
%         VSs = PhaseData(:,ir,2);
%         RSs = PhaseData(:,ir,3);
        
        %     subplot(2,3,ir);
        %     polarplot(deg2rad(Phs),RSs./(VSs.^2)./2,'ok')
        %     set(gca,'Color','none','RAxisLocation',130)
        %     set(gca,'ThetaZeroLocation','bottom',...
        %         'ThetaDir','clockwise');
        %     set(gca,'rlim',[0 8000])
        
        histplot_phs = zeros(1,numel(radvector));
        histplot_sum = zeros(1,numel(radvector));
        histplot_avg = zeros(1,numel(radvector));
        for irad = 1:numel(radvector)
            if irad==1 || irad==numel(radvector)
                theseUns =  Phs>radvector(end-1) | Phs<=radvector(2);
            else
                theseUns =  Phs>radvector(irad-1) & Phs<=radvector(irad+1);
            end
            
            histplot_phs(irad) = sum(theseUns);
            if sum(theseUns)>0
                histplot_sum(irad) = sum(RSs(theseUns)./(VSs(theseUns).^2)./2);
                histplot_avg(irad) = sum(RSs(theseUns)./(VSs(theseUns).^2)./2)/sum(theseUns);
            end
        end
        
        figure(hfsum); hold on
        
        subplot(3,5,ir);
        polarplot(deg2rad(radvector),histplot_phs,'Color',[0.173 0.171 0.686],'LineWidth',2)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        rmax = 36;
        set(gca,'rlim',[0 rmax],'rtick',linspace(0,rmax,4),'rticklabel',{'' '' '' num2str(rmax)})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        if ir==1
            title('Phase distribution, sync units')
        else
            title([num2str(AMrates(ir)) ' Hz'])
        end
        
        
        % Average number of spikes per unit
        subplot(3,5,ir+5);
        polarplot(deg2rad(radvector),histplot_avg,'Color',[220 0 0]./255,'LineWidth',2)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        rmax = 1000;
        set(gca,'rlim',[0 rmax],'rtick',linspace(0,rmax,4),'rticklabel',{'' '' '' num2str(rmax)})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        if ir==1
            title('Avg N spks per unit')
        end
        
        subplot(3,5,ir+10);
        polarplot(deg2rad(radvector),histplot_sum,'Color',[150 0 200]./255,'LineWidth',2)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        rmax = 8000;
        set(gca,'rlim',[0 rmax],'rtick',linspace(0,rmax,4),'rticklabel',{'' '' '' num2str(rmax)})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        if ir==1
            title('Total N spks, sync units')
        end
        
        
    end
    
    % Save figure
    figure(hfsum); hold on
    set(hfsum,'PaperOrientation','landscape')
    print(fullfile(savedir,'Phs_Nspk_Distributions_sync'),'-dpdf','-bestfit')
    
end



% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% Look at real time MPH data
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load(fullfile(fn.figs,'PopMPH','MPHdata.mat'))

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end


% Data settings
switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
%         Boundaries  = [-1 round(10.^[0.5]) ];
end


%% Max Peak: time x duration

if plotF3
    
    for iC = CellTypes
        
        % Data settings
        switch iC{:}
            case 'RS'
                theseCells = iRS;
            case 'NS'
                theseCells = iNS;
        end
        
        switch useFR
            case 'z'
                ThisData    = zFR_vec(theseCells,:,:);
                Boundaries  = [-1 0 0.25 0.5 1 2];
            case 'log'
                ThisData    = FR_vec(theseCells,:,:);
                Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
                %         Boundaries  = [-1 round(10.^[0.5]) ];
        end
        
        
        hfPh=figure;
        set(hfPh,'Position',fullscreen)
        
        t_onset   = nan(size(ThisData,1),3);
        t_offset  = nan(size(ThisData,1),3);
        t_up      = nan(size(ThisData,1),3);
        t_mid     = nan(size(ThisData,1),3);
        t_Lat     = nan(size(ThisData,1),3);
        pk_height = nan(size(ThisData,1),3);
        pk_time   = nan(size(ThisData,1),3);
        
        for iUn = 1:numel(theseCells)
            
            if sum(UnitData(theseCells(iUn)).VSdata_spk(2,2:4)>13.1)<3
                continue
            end
            
            for ist = 1:3
                
                %         plot(ThisData(iUn,:,ist),'Color',colors(ist,:),'LineWidth',3)
                %         hold on
                
                pdms = ceil(1000/AMrates(ist));
                
                halfMax = ( max(ThisData(iUn,:,ist)) - min(ThisData(iUn,:,ist)) ) /2;
                [PKS,LOCS] = findpeaks(ThisData(iUn,:,ist),'MinPeakHeight',halfMax);
                
                if isempty(PKS)
                    continue
                end
                
                % ASSUME just one peak
                [~,ipk] = max(PKS);
                pk_height(iUn,ist) = PKS(ipk);
                pk_time(iUn,ist) = LOCS(ipk);
                
                % Find time before surpassing halfMax, and time after falling below
                for ims = [LOCS(ipk):-1:1 pdms:-1:LOCS(ipk)]
                    if ThisData(iUn,ims,ist)<halfMax
                        t_onset(iUn,ist) = ims;
                        break
                    end
                end
                
                for ims = [LOCS(ipk):pdms 1:LOCS(ipk)]
                    if ThisData(iUn,ims,ist)<halfMax
                        t_offset(iUn,ist) = ims;
                        break
                    end
                end
                
                % Calculate peak duration
                if t_offset(iUn,ist)>=t_onset(iUn,ist)
                    t_up(iUn,ist) = t_offset(iUn,ist)-t_onset(iUn,ist);
                else
                    t_up(iUn,ist) = t_offset(iUn,ist)+pdms-t_onset(iUn,ist);
                end
                
                % Time of midpoint of peak
                t_mid(iUn,ist) = t_onset(iUn,ist)+t_up(iUn,ist)/2;
                if t_mid(iUn,ist)>pdms
                    t_mid(iUn,ist) = t_mid(iUn,ist)-pdms;
                end
                
                % Set latency value
                switch whichLat 
                    case 'Onset'
                        t_Lat(iUn,ist) = t_onset(iUn,ist);
                    case 'Mid'
                        t_Lat(iUn,ist) = t_mid(iUn,ist);
                    case 'Offset'
                        t_Lat(iUn,ist) = t_offset(iUn,ist);
                end
                
                % Add to plot
                subplot(2,3,ist); hold on
                plot(t_Lat(iUn,ist),t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
                
            end %ist
        end %iUn
        
        % Finish plots
        for ist = 1:3
            pdms = ceil(1000/AMrates(ist));
            
            subplot(2,3,ist); hold on
            axis square
            xlim([0 pdms])
            ylim([0 pdms])
            set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
            xlabel([whichLat ' time (ms)'])
            ylabel('Peak duration (ms)')
            title([num2str(AMrates(ist)) ' Hz'])
        end
        
        
        % Below, plot relationship between metrics of consecutive rates
        for ist = 1:2
            pdms = ceil(1000/AMrates(ist));
            
            subplot(2,3,3+ist); cla
            plot([0 pdms],[0 pdms/2],'Color',0.5*[1 1 1])
            hold on
            plot(t_Lat(:,ist),t_Lat(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
            
            axis square
            
            xlim([0 pdms])
            ylim([0 pdms])
            set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
            xlabel([whichLat ' time, ' num2str(AMrates(ist)) ' Hz'])
            ylabel([whichLat ' time, ' num2str(AMrates(ist+1)) ' Hz'])
        end
        
        % Durations
        alt_cols = [colors(1,[2 1 3]); colors(2,[3 2 1])];
        
        for ist = 1:2
            pdms = ceil(1000/AMrates(ist));
            
            subplot(2,3,6);
            plot([0 pdms],[0 pdms/2],'Color',alt_cols(ist,:))
            hold on
            plot(t_up(:,ist),t_up(:,ist+1),'.','Color',alt_cols(ist,:),'MarkerSize',20)
            
        end
        
        pdms = ceil(1000/AMrates(1));
        plot([0 pdms],[0 pdms],'Color',0.5*[1 1 1])
        axis square
        xlim([0 pdms])
        ylim([0 pdms])
        set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
        xlabel(['Peak duration, ' num2str(AMrates(ist)) ' Hz'])
        ylabel(['Peak duration, ' num2str(AMrates(ist+1)) ' Hz'])
        
        suptitle(iC{:})
        
        print_eps_kp(gcf,fullfile(savedir,['Peak_' whichLat '_' iC{:}]))
        
    end %CellTypes
end


%%

if plotF5
    
    hfPh=figure;
    set(hfPh,'Position',fullscreen)
    
    for iC = CellTypes
        
        % Data settings
        switch iC{:}
            case 'RS'
                theseCells = iRS;
                ispr = 1;
            case 'NS'
                theseCells = iNS;
                ispr = 2;
        end
        
        
        % Get phases
        allPhases = nan(numel(theseCells),3);
        for iUn = 1:numel(theseCells)
            allPhases(iUn,:) = UnitData(theseCells(iUn)).Phase_spk(2:4);
        end
        
        
        switch useFR
            case 'z'
                ThisData    = zFR_vec(theseCells,:,:);
                Boundaries  = [-1 0 0.25 0.5 1 2];
            case 'log'
                ThisData    = FR_vec(theseCells,:,:);
                Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
                %         Boundaries  = [-1 round(10.^[0.5]) ];
        end
        
        t_onset   = nan(size(ThisData,1),3);
        t_offset  = nan(size(ThisData,1),3);
        t_up      = nan(size(ThisData,1),3);
        t_mid     = nan(size(ThisData,1),3);
        pk_height = nan(size(ThisData,1),3);
        pk_time   = nan(size(ThisData,1),3);
        
        for iUn = 1:numel(theseCells)
            
            if sum(UnitData(theseCells(iUn)).VSdata_spk(2,2:4)>13.1)<3
                continue
            end
            
            for ist = 1:3
                
                %         plot(ThisData(iUn,:,ist),'Color',colors(ist,:),'LineWidth',3)
                %         hold on
                
                pdms = ceil(1000/AMrates(ist));
                
                halfMax = ( max(ThisData(iUn,:,ist)) - min(ThisData(iUn,:,ist)) ) /2;
                [PKS,LOCS] = findpeaks(ThisData(iUn,:,ist),'MinPeakHeight',halfMax);
                
                if isempty(PKS)
                    continue
                end
                
                % ASSUME just one peak
                [~,ipk] = max(PKS);
                pk_height(iUn,ist) = PKS(ipk);
                pk_time(iUn,ist) = LOCS(ipk);
                
                % Find time before surpassing halfMax, and time after falling below
                for ims = [LOCS(ipk):-1:1 pdms:-1:LOCS(ipk)]
                    if ThisData(iUn,ims,ist)<halfMax
                        t_onset(iUn,ist) = ims;
                        break
                    end
                end
                
                for ims = [LOCS(ipk):pdms 1:LOCS(ipk)]
                    if ThisData(iUn,ims,ist)<halfMax
                        t_offset(iUn,ist) = ims;
                        break
                    end
                end
                
                % Calculate peak duration
                if t_offset(iUn,ist)>=t_onset(iUn,ist)
                    t_up(iUn,ist) = t_offset(iUn,ist)-t_onset(iUn,ist);
                else
                    t_up(iUn,ist) = t_offset(iUn,ist)+pdms-t_onset(iUn,ist);
                end
                
                % Time of midpoint of peak
                t_mid(iUn,ist) = t_onset(iUn,ist)+t_up(iUn,ist)/2;
                if t_mid(iUn,ist)>pdms
                    t_mid(iUn,ist) = t_mid(iUn,ist)-pdms;
                end
                
                
                % Add to plot
                subplot(2,3,ist+(ispr-1)*3); hold on
                plot(t_mid(iUn,ist),allPhases(iUn,ist)/360*pdms,'.','Color',colors(ist,:),'MarkerSize',10)
                
            end %ist
        end %iUn
        
        % Finish plots
        for ist = 1:3
            pdms = ceil(1000/AMrates(ist));
            
            subplot(2,3,ist+(ispr-1)*3); hold on
            plot([0 pdms],[0 pdms],'-','Color',0.6.*[1 1 1])
            
            axis square
            xlim([0 pdms])
            ylim([0 pdms])
            set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
            xlabel('Mid time (ms)')
            ylabel('Mean phase (ms)')
            title([num2str(AMrates(ist)) ' Hz'])
        end
        
    end %CellTypes
    
    suptitle(CellTypes)
    
    print_eps_kp(gcf,fullfile(savedir,'MidTime_Phase'))

end


end



