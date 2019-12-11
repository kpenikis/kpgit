function MoreRespChar
% 
% MoreRespChar
%
%  Group cells in population by their MPHs. 
% 

global AMrates useFR Boundaries

%~~~~~~~~~~~~~~~~~~~~~
plotF1   = 0;
plotF2   = 0;
plotF3   = 0;
plotF4   = 1;
%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
useGrp   =   'peakFR'; 'dynRange'; 'tuning'; 'phase'; 
%~~~~~~~~~~~~~~~~~~~~~
ngrps    = 5;
AMrates = [2 4 8 16 32];

% close all

%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


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
       
savedir = fullfile(fn.figs,'RespChar');
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
    
hfPk=figure; 
set(hfPk,'Position',fullscreen)

% show non sync units?
% plot cumulative nspikes behind? 

t_onset   = nan(size(ThisData,1),3);
t_offset  = nan(size(ThisData,1),3);
t_up      = nan(size(ThisData,1),3);
t_mid     = nan(size(ThisData,1),3);
pk_height = nan(size(ThisData,1),3);
pk_time   = nan(size(ThisData,1),3);

for iUn = 1:size(ThisData,1)
    
    if sum(UnitData(iUn).VSdata_spk(2,2:4)>13.1)<3
        continue
    end
    
%     pause(3)
%     clf
    
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
        subplot(2,3,ist); hold on
%         plot(t_onset(iUn,ist), t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
        plot(t_mid(iUn,ist),t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
        
%         plot(t_onset(iUn,ist)+[0 t_up(iUn,ist)],[halfMax halfMax],'Color',colors(ist,:),'LineWidth',6)
        
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
    xlabel('Mid time (ms)')
    ylabel('Peak duration (ms)')
    title([num2str(AMrates(ist)) ' Hz'])
end


% Below, plot relationship between metrics of consecutive rates
for ist = 1:2
    pdms = ceil(1000/AMrates(ist));
    
    subplot(2,3,3+ist); cla
    plot([0 pdms],[0 pdms/2],'Color',0.5*[1 1 1])
    hold on
    plot(t_mid(:,ist),t_mid(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
%     plot(t_onset(:,ist),t_onset(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
    
    axis square
    
    xlim([0 pdms])
    ylim([0 pdms])
    set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
    xlabel(['Mid time, ' num2str(AMrates(ist)) ' Hz'])
    ylabel(['Mid time, ' num2str(AMrates(ist+1)) ' Hz'])
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


savedir = fullfile(fn.figs,'RespChar');
print_eps_kp(gcf,fullfile(savedir,'Peak_Mids'))

end




%% Time x duration, latency from peak time

if plotF4 
    
Rotate =0;
Reverse=1;

hfPk=figure; 
set(hfPk,'Position',fullscreen)

% show non sync units?
% plot cumulative nspikes behind? 

t_onset   = nan(size(ThisData,1),3);
t_offset  = nan(size(ThisData,1),3);
t_up      = nan(size(ThisData,1),3);
t_mid     = nan(size(ThisData,1),3);
pk_height = nan(size(ThisData,1),3);
pk_time   = nan(size(ThisData,1),3);

for iUn = 1:size(ThisData,1)
    
    if sum(UnitData(iUn).VSdata_spk(2,2:4)>13.1)<3
        continue
    end
    
%     pause(3)
%     clf
    
    for ist = 1:3
        
        Data=[];
        PKS =[];
        LOCS=[];
        
        pdms = 1000/AMrates(ist);
        
        if Rotate
            Data = [ ThisData(iUn,floor(pdms/2):floor(pdms),ist) ThisData(iUn,1:floor(pdms/2)-1,ist)];
        elseif Reverse
            Data = fliplr(ThisData(iUn,1:floor(pdms),ist));
        end
%         plot(ThisData(iUn,:,ist),'Color',colors(ist,:),'LineWidth',3)
%         hold on
        
        pdms = floor(1000/AMrates(ist));
        
        halfMax = ( max(Data) - min(Data) ) /2;
        [PKS,LOCS] = findpeaks(Data,'MinPeakHeight',halfMax);
        
        if isempty(PKS)
            continue
        end
        
        % ASSUME just one peak
        [~,ipk] = max(PKS);
        pk_height(iUn,ist) = PKS(ipk);
        pk_time(iUn,ist) = LOCS(ipk);
        
        % Find time before surpassing halfMax, and time after falling below        
        for ims = [LOCS(ipk):-1:1 pdms:-1:LOCS(ipk)]
            if Data(ims)<halfMax
                t_onset(iUn,ist) = ims;
                break
            end
        end
        if t_onset(iUn,ist)>pdms
            keyboard
        end
        
        for ims = [LOCS(ipk):pdms 1:LOCS(ipk)]
            if Data(ims)<halfMax
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
        subplot(2,3,ist); hold on
%         plot(t_onset(iUn,ist), t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
%         plot(t_mid(iUn,ist),t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
        plot(t_offset(iUn,ist),t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
        
%         plot(t_onset(iUn,ist)+[0 t_up(iUn,ist)],[halfMax halfMax],'Color',colors(ist,:),'LineWidth',6)
        
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
    xlabel('Offset time (ms)')
    ylabel('Peak duration (ms)')
    title([num2str(AMrates(ist)) ' Hz'])
end


% Below, plot relationship between metrics of consecutive rates
for ist = 1:2
    pdms = floor(1000/AMrates(ist));
    
    subplot(2,3,3+ist); cla
    plot([0 pdms],[0 pdms/2],'Color',0.5*[1 1 1])
    hold on
    plot(t_offset(:,ist),t_offset(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
%     plot(t_mid(:,ist),t_mid(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
%     plot(t_onset(:,ist),t_onset(:,ist+1),'.','Color',0.3*[1 1 1],'MarkerSize',20)
    
    axis square
    
    xlim([0 pdms])
    ylim([0 pdms])
    set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
    xlabel(['Offset time, ' num2str(AMrates(ist)) ' Hz'])
    ylabel(['Offset time, ' num2str(AMrates(ist+1)) ' Hz'])
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


savedir = fullfile(fn.figs,'RespChar');
print_eps_kp(gcf,fullfile(savedir,'Peak_Offsets_reversed'))


end







keyboard







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
set(gcf,'Position',sz1)
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

