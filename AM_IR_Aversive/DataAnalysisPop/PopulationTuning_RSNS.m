function PopulationTuning_RSNS
% Figure 3

% [ ] proportion units significantly above and below silence baseline 
% [ ] distribution of BMF-VS
% [ ] cdf to pdf, don't plot z-score
% [ ] trial variability plots 
% [ ] bar: mean FF as function of stimulus (last 750 ms)
% [ ] scatter: FF as a function of FR  (last 750 ms)
% [ ] median CorrScore as a function of stimulus (last 750 ms)
% [ ] mean percentage ?miss? trials as a function of stimulus (last 750 ms)


alfaVS = 0.0002;
alfaFR = 0.01;

minTrs  = 12;

AMrates = [2 4 8 16 32];
Stimuli = [9 1:8];

binsmth = 10;
getMPH  = 0;

% Set colors
colors = [ 150 150 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        
% Load Unit data files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

[~,UnitData] = identifyResponsiveUnits(UnitData);

%%

VS_sig   = zeros(numel(UnitData),numel(AMrates));
BMF_VS   = zeros(numel(UnitData),1);

FR_sig05 = zeros(numel(UnitData),numel(AMrates));
FR_sig01 = zeros(numel(UnitData),numel(AMrates));
FR_dp    = zeros(numel(UnitData),numel(AMrates));
BMF_FR   = zeros(numel(UnitData),1);

FRmins   = nan(numel(UnitData),numel(AMrates));
FRmaxs   = nan(numel(UnitData),numel(AMrates));

FFmeans  = nan(numel(UnitData),numel(AMrates)+1);

zFRs     = nan(numel(UnitData),7);
dFRs     = nan(numel(UnitData),8);

PdcRespType = cell(numel(UnitData),5);

for iUn = 1:numel(UnitData)
    
    % Percent of units synchronized by AMrate
    VS_sig(iUn,:) = UnitData(iUn).VSdata_spk(3,2:6)<alfaVS;
    
    % Distribution of BMF-VS and N ns
    if ~isempty(UnitData(iUn).iBMF_VS)
        BMF_VS(iUn) = UnitData(iUn).iBMF_VS;
    end
    
    % Percent of units FR > silence per stimulus, dprime
    FR_dp(iUn,:) = abs(UnitData(iUn).dp_mat(2:6,2))';
    
    % Distribution of BMF-FR and N ns
    if ~isempty(UnitData(iUn).iBMF_FR)
        BMF_FR(iUn) = UnitData(iUn).iBMF_FR;
    end    
    
    % Avg FR per stimulus, z-scored
    zFRs(iUn,:) = UnitData(iUn).FR_nrm(2:end);
    
    
%     % Filter FRtrials to datapoints above minTrs
%     if any(sum(~isnan(UnitData(iUn).FR_raw_tr(:,2:end)),1)<minTrs & sum(~isnan(UnitData(iUn).FR_raw_tr(:,2:end)),1)>1)
%         keyboard
%     end
%     dFRs(iUn,:) = mean(UnitData(iUn).FR_raw_tr(:,2:end),1,'omitnan') - UnitData(iUn).BaseFR;
%     
    
    % - - - -  The following calculations require raw spike data - - - - - 
    clear FRtrials
    get_trial_data_posthoc
    
    
    % Filter FRtrials to datapoints above minTrs
    if any( sum(~isnan(FRtrials),1)<minTrs & sum(~isnan(FRtrials),1)>1 )
        idx=find( sum(~isnan(FRtrials),1)<minTrs & sum(~isnan(FRtrials),1)>1 );
        FRtrials(:,idx) = nan;
        fprintf('removed stim %i \n',idx)
    end
    
    mFRs(iUn,:) = mean(FRtrials(:,2:end),1,'omitnan');
    
    % New way to get dFRs
    dFRs(iUn,:) = mean(FRtrials(:,2:end),1,'omitnan') - UnitData(iUn).BaseFR;
    
    
    % Percent of units FR significantly > silence
    pvals = nan(1,5);
    for ist = 1:5
        try
            pvals(ist) = ranksum(FRtrials(:,1),FRtrials(:,ist+2));
        end
    end
    
    FRsign = sign(median(FRtrials(:,3:7),'omitnan')-median(FRtrials(:,1),'omitnan'));
%     FRsign = sign(mean(FRtrials(:,3:7),'omitnan')-mean(FRtrials(:,1),'omitnan'));
    
    [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,0.05);
    FR_sig05(iUn,:) = iIR_sig_bonferroni.*FRsign;
    
    [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,0.01);
    FR_sig01(iUn,:) = iIR_sig_bonferroni.*FRsign;
    
    
    % Fano factor for each stimulus (last 750 ms)
%     FFmeans(iUn,:) = var(FRtrials,1)./mean(FRtrials,1);
    
    % Range of FR across trials
%     FRmins(iUn,:) = min(FRtrials(:,2:end),[],1); %min(UnitData(iUn).FR_raw_tr(:,2:6))
%     FRmaxs(iUn,:) = max(FRtrials(:,2:end),[],1); %max(UnitData(iUn).FR_raw_tr(:,2:6))
    
    
    % Get Pdc response type                      [ early  late  p_val ]
    for irate = 1:5
        PdcRespType{iUn,irate} = 'n';
        if ~isempty(UnitData(iUn).DeltaNspk) && ~isempty(UnitData(iUn).DeltaNspk{irate})
            if (UnitData(iUn).DeltaNspk{irate}(3))<alfaFR && (UnitData(iUn).DeltaNspk{irate}(1) > UnitData(iUn).DeltaNspk{irate}(2))
                % adapting
                PdcRespType{iUn,irate} = 'A';
            elseif (UnitData(iUn).DeltaNspk{irate}(3))<alfaFR && (UnitData(iUn).DeltaNspk{irate}(1) < UnitData(iUn).DeltaNspk{irate}(2))
                % facilitating
                PdcRespType{iUn,irate} = 'F';
            else
                PdcRespType{iUn,irate} = 'S';
            end
        end
    end
end %iUn

BMF_FR(BMF_FR==0) = nan;
BMF_VS(BMF_VS==0) = nan;


%% Plot results

set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallrect    = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


% CellTypes
iRS = UnitInfo.TroughPeak>0.43;
iNS = UnitInfo.TroughPeak<=0.43;


%% Adapted/facilitated

fprintf('\n RS cells: \n')
fprintf('A\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir),sum(strcmp(PdcRespType(iRS,ir),'A'))/sum(iRS)*100)
end
fprintf('\nF\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir),sum(strcmp(PdcRespType(iRS,ir),'F'))/sum(iRS)*100)
end
fprintf('\nall\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir), ...
        (sum(strcmp(PdcRespType(iRS,ir),'A'))+sum(strcmp(PdcRespType(iRS,ir),'F'))) /sum(iRS)*100 )
end

fprintf('\n \n NS cells: \n')
fprintf('A\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir),sum(strcmp(PdcRespType(iNS,ir),'A'))/sum(iNS)*100)
end
fprintf('\nF\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir),sum(strcmp(PdcRespType(iNS,ir),'F'))/sum(iNS)*100)
end
fprintf('\nall\n')
for ir = 1:5
    fprintf(' %iHz: %0.1f  ',AMrates(ir), ...
        (sum(strcmp(PdcRespType(iNS,ir),'A'))+sum(strcmp(PdcRespType(iNS,ir),'F'))) /sum(iNS)*100 )
end

fprintf('\n')


%% PLOT

ycountmax = 60; %10*ceil(sum(iRS)/20);

hf = figure;
set(hf,'Position',tallrect)

%--------------------------------------------------------------------------
% FR probability distributions
%--------------------------------------------------------------------------

FRmin = .001;
FRmax = 1000;
FRstp =   25;

%~~~~~~~~~~~~~~ RS cells
subplot(5,2,1)
%   [X] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
for i=1
    
% First, distribution for Silence
data = round([UnitData(iRS).BaseFR]',2);
[values,edges] = histcounts(log10(data),linspace(log10(FRmin),log10(FRmax),FRstp));
x = (edges(1:end-1) + mode(diff(edges))/2);
f = fit(x',values','gauss1');
xnew = linspace(log10(FRmin),log10(FRmax),100);
z = feval(f,xnew);
plot(xnew,z./sum(z),'k','LineWidth',2)
% plot(x,values./sum(values),'k','LineWidth',2)

% Save distribution for stats
FRdists_fit = nan(length(z),1+size(dFRs,2));
FRdists_dat = nan(length(data),1+size(dFRs,2));
FRdists_fit(:,1) = z./sum(z);
FRdists_dat(:,1) = data;

DataVec_FR = data;
DataVec_st = ones(size(data));
DataVec_cl = ones(size(data));

hold on
for ist = 1:size(dFRs,2)
    
    data = round( dFRs((iRS),ist)+[UnitData(iRS).BaseFR]' ,2);
    
    [values,edges] = histcounts(log10(data),linspace(log10(FRmin),log10(FRmax),FRstp));
    x = (edges(1:end-1) + mode(diff(edges))/2);
    
    f = fit(x',values','gauss1');
    xnew = linspace(log10(FRmin),log10(FRmax),100);
    z=[]; z = feval(f,xnew);
    plot(xnew,z./sum(z),'-','Color',colors(ist,:),'LineWidth',2);
%     plot(x,values./sum(values),'-','Color',colors(ist,:),'LineWidth',2);
    
    FRdists_fit(:,ist+1) = z./sum(z);
    FRdists_dat(:,ist+1) = data;
    
    DataVec_FR = [DataVec_FR; data];
    DataVec_st = [DataVec_st; (ist+1)*ones(size(data))];
    DataVec_cl = [DataVec_cl; ones(size(data))];
    
end

plot(log10(median(FRdists_dat(:),1,'omitnan')),0,'o','MarkerFaceColor','k','MarkerEdgeColor','none')

% axis square
set(gca,'Color','none','ytick',[0 0.12])
ylim([0 0.12])
xlabel('log(FR) (sp/s)')
ylabel('Probability')
grid on
title('RS')

% p_dist = anova1(FRdists_fit);
% [p_data,~,statss] = anova1(FRdists_dat);
% multcompare(statss);
% title(sprintf('RS, F=%0.2f',p_data))
end


%~~~~~~~~~~~~~~ NS cells
figure(hf); hold on
subplot(5,2,2)
%   [ ] [X]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
for i=1
    
% First, distribution for Silence
data = round([UnitData(iNS).BaseFR]',2);
[values,edges] = histcounts(log10(data),linspace(log10(0.001),log10(1000),FRstp));
x = (edges(1:end-1) + mode(diff(edges))/2);
f = fit(x',values','gauss1');
xnew = linspace(log10(0.001),log10(1000),100);
z = feval(f,xnew);
plot(xnew,z./sum(z),'k','LineWidth',2)
% plot(x,values./sum(values),'k','LineWidth',2)

% Save distribution for stats
FRdists_fit = nan(length(z),1+size(dFRs,2));
FRdists_dat = nan(length(data),1+size(dFRs,2));
FRdists_fit(:,1) = z./sum(z);
FRdists_dat(:,1) = data;

DataVec_FR = [DataVec_FR; data];
DataVec_st = [DataVec_st; ones(size(data))];
DataVec_cl = [DataVec_cl; 2*ones(size(data))];

hold on
for ist = 1:size(dFRs,2)
    
    data = round( dFRs((iNS),ist)+[UnitData(iNS).BaseFR]' ,2);
    
    [values,edges] = histcounts(log10(data),linspace(log10(0.001),log10(1000),FRstp));
    x = (edges(1:end-1) + mode(diff(edges))/2);
    
    f = fit(x',values','gauss1');
    xnew = linspace(log10(0.001),log10(1000),100);
    z=[]; z = feval(f,xnew);
    plot(xnew,z./sum(z),'-','Color',colors(ist,:),'LineWidth',2);
%     plot(x,values./sum(values),'-','Color',colors(ist,:),'LineWidth',2);
%     plot(x,medfilt1(values./sum(values)),'-','Color',colors(ist,:),'LineWidth',2);
    
    FRdists_fit(:,ist+1) = z./sum(z);
    FRdists_dat(:,ist+1) = data;
    
    DataVec_FR = [DataVec_FR; data];
    DataVec_st = [DataVec_st; (ist+1)*ones(size(data))];
    DataVec_cl = [DataVec_cl; 2*ones(size(data))];
end

plot(log10(median(FRdists_dat(:),1,'omitnan')),0,'o','MarkerFaceColor','k','MarkerEdgeColor','none')

% axis square
set(gca,'Color','none','ytick',[0 0.12])
ylim([0 0.12])
% xlabel('log(FR) (sp/s)')
% ylabel('Probability')
grid on
title('NS')

end

% p_dist = anova1(FRdists_fit);
% [p_data,~,statss] = anova1(FRdists_dat);
% multcompare(statss);
% title(sprintf('NS, F=%0.2f',p_data))


% [p_data,~,statss] = anovan(log10(DataVec_FR+1e-6),[DataVec_st DataVec_cl]);
% [p_data,~,statss] = kruskalwallis(log10(DataVec_FR),DataVec_st);


%%% One cell type at a time

% RS:   16 and 32 Hz < DB

[p_kw_RS,~,stats_RS] = kruskalwallis(log10(DataVec_FR(DataVec_cl==1)),DataVec_st(DataVec_cl==1),'off');
figure;
multcompare(stats_RS,'Dimension',1,'CType','bonferroni');
% title('hsd') %bonferroni
fprintf('  RS: kw p=%0.4f\n',p_kw_RS)


% NS:   NONE  or  32 Hz < 8 Hz

[p_kw_NS,~,stats_NS] = kruskalwallis(DataVec_FR(DataVec_cl==2),DataVec_st(DataVec_cl==2),'off');
figure;
multcompare(stats_NS,'Dimension',1,'CType','bonferroni');
% title('hsd') %scheffe
fprintf('  NS: kw p=%0.4f\n',p_kw_NS)



%% Rate response analyses

figure(hf); hold on

%--------------------------------------------------------------------------
%                     Percent different than spont.
%--------------------------------------------------------------------------
%~~~~~~~~~~~~~~ RS cells
subplot(5,2,3)
%   [ ] [ ]
%   [X] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
for i=1

% ranksum test
plot(1:5,sum(FR_sig05(iRS,:)>0,1)/size(FR_sig05(iRS,:),1),'-','LineWidth',4,'Color',[0.3 0.5 0.95])
hold on
plot(1:5,sum(FR_sig05(iRS,:)~=0,1)/size(FR_sig05(iRS,:),1),'-','LineWidth',4,'Color',[0.2 0.33 0.8])
% plot(1:5,sum(FR_sig01(iRS,:),1)/size(FR_sig01(iRS,:),1),'-','LineWidth',4,'Color',[0.2 0.33 0.8])
set(gca,'Color','none','xtick',1:5,'xticklabel',AMrates,'tickdir','out',...
    'ticklength',[0.02 0.02],'ytick',[0 1])
xlim([0 6])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Prop. cells')
% axis square
box on
end

%~~~~~~~~~~~~~~ NS cells
subplot(5,2,4)
%   [ ] [ ]
%   [ ] [X]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
for i=1
    
% ranksum test
plot(1:5,sum(FR_sig05(iNS,:)>0,1)/size(FR_sig05(iNS,:),1),'-','LineWidth',4,'Color',[0.3 0.5 0.95])
hold on
plot(1:5,sum(FR_sig05(iNS,:)~=0,1)/size(FR_sig05(iNS,:),1),'-','LineWidth',4,'Color',[0.2 0.33 0.8])
% plot(1:5,sum(FR_sig01(iNS,:),1)/size(FR_sig01(iNS,:),1),'-','LineWidth',4,'Color',[0.2 0.33 0.8])
set(gca,'Color','none','xtick',1:5,'xticklabel',AMrates,'tickdir','out',...
    'ticklength',[0.02 0.02],'ytick',[0 1])
xlim([0 6])
ylim([0 1])
xlabel('AM rate (Hz)')
% ylabel('Proportion FR>silence')
% axis square
box on
end

%--------------------------------------------------------------------------
%                          BMF distributions
%--------------------------------------------------------------------------

%~~~~~~~~~~~~~~ RS cells
subplot(5,2,5)
%   [ ] [ ]
%   [ ] [ ]
%   [X] [ ]
%   [ ] [ ]
%   [ ] [ ]
for i=1
    
histogram(BMF_FR(iRS),'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax])
xlim([0 6])
ylim([0 ycountmax])
xlabel('BMF (FR)')
ylabel('# units')
% axis square
box on
end

%~~~~~~~~~~~~~~ NS cells
subplot(5,2,6)
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [X]
%   [ ] [ ]
%   [ ] [ ]
for i=1
histogram(BMF_FR(iNS),'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax])
xlim([0 6])
ylim([0 ycountmax/2])
xlabel('BMF (FR)')
% ylabel('NS # units')
% axis square
box on
end


%% Vector strength analyses

%--------------------------------------------------------------------------
%                       Percent synchronized
%--------------------------------------------------------------------------

%~~~~~~~~~~~~~~ RS cells
subplot(5,2,7)
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [X] [ ]
%   [ ] [ ]

plot(1:5,sum(VS_sig(iRS,:),1)/size(VS_sig(iRS,:),1),'.-r','MarkerSize',30,'LineWidth',4)
set(gca,'Color','none','xtick',1:5,'xticklabel',AMrates,'tickdir','out',...
    'ticklength',[0.02 0.02],'ytick',[0 1])
xlim([0 6])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Proportion sig sync')
% axis square
box on

%~~~~~~~~~~~~~~ NS cells
subplot(5,2,8)
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [X]
%   [ ] [ ]
plot(1:5,sum(VS_sig(iNS,:),1)/size(VS_sig(iNS,:),1),'.-r','MarkerSize',30,'LineWidth',4)
set(gca,'Color','none','xtick',1:5,'xticklabel',AMrates,'tickdir','out',...
    'ticklength',[0.02 0.02],'ytick',[0 1])
xlim([0 6])
ylim([0 1])
xlabel('AM rate (Hz)')
% ylabel('Proportion sig sync')
% axis square
box on


%--------------------------------------------------------------------------
%                          BMF distributions
%--------------------------------------------------------------------------

%~~~~~~~~~~~~~~ RS cells
subplot(5,2,9)
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [X] [ ]

histogram(BMF_VS(iRS),'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax])
xlim([0 6])
ylim([0 ycountmax])
xlabel('BMF (VS)')
ylabel('# units')
% axis square
box on

%~~~~~~~~~~~~~~ RS cells
subplot(5,2,10)
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [ ]
%   [ ] [X]

histogram(BMF_VS(iNS),'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax])
xlim([0 6])
ylim([0 ycountmax/2])
xlabel('BMF (VS)')
% ylabel('N units')
% axis square
box on



% Save figure
savedir = fullfile(fn.figs,'PopulationTuning','Apr2020');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,'TuningSummary_RSNS_fitLogNorm'))


%% Overlap of FR and Sync responses

figure;
set(gcf,'Position',fullscreen)

fileID = fopen(fullfile(savedir,'Overlap_SyncFR_RSNS.txt'),'w');

% RS cells
fprintf(fileID,'\n----- RS cells -----\n');
for ir = 1:5
    
    % Sync
    r1_S1_FRnc = FR_sig05(iRS,ir)==0  & VS_sig(iRS,ir)==1;
    r3_S1_FR0  = FR_sig05(iRS,ir)==-1 & VS_sig(iRS,ir)==1;
    r5_S1_FR1  = FR_sig05(iRS,ir)==1  & VS_sig(iRS,ir)==1;
    
    % Non-sync
    r2_S0_FRnc = FR_sig05(iRS,ir)==0  & VS_sig(iRS,ir)==0;
    r4_S0_FR0  = FR_sig05(iRS,ir)==-1 & VS_sig(iRS,ir)==0;
    r6_S0_FR1  = FR_sig05(iRS,ir)==1  & VS_sig(iRS,ir)==0;
    
    nsync  = sum(VS_sig(iRS,ir)==1);
    nNsync = sum(VS_sig(iRS,ir)==0);
    vertmid = round(100*sum(VS_sig(iRS,ir)==1)/sum(iRS),1);
    
    subplot(2,5,ir)
    xlim([0 100])
    ylim([0 100])
    hold on
    axis square
    set(gca,'Color','none')
    box off
    
    % Sync
    r1_x = [0 vertmid vertmid 0];
    r1_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1) 100 100];
    fill(r1_x,r1_y,0.6.*[1 1 1],'EdgeColor','none')
    
    r3_x = [0 vertmid vertmid 0];
    r3_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)...
        100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1)];
    fill(r3_x,r3_y,[0.9 0.2 0.2],'EdgeColor','none')
    
    r5_x = [0 vertmid vertmid 0];
    r5_y = [0 0 ...
        100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)];
    fill(r5_x,r5_y,[0.2 0.2 0.9],'EdgeColor','none')
    
    % Non-sync
    r2_x = [vertmid 100 100 vertmid];
    r2_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100 100];
    fill(r2_x,r2_y,0*[1 1 1],'EdgeColor','none')
    
    r4_x = [vertmid 100 100 vertmid];
    r4_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1)];
    fill(r4_x,r4_y,[0.6 0 0],'EdgeColor','none')
    
    r6_x = [vertmid 100 100 vertmid];
    r6_y = [0 0 ...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)];
    fill(r6_x,r6_y,[0 0 0.6],'EdgeColor','none')
    
    
    % Print results
    fprintf(fileID,'%i Hz: \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n',...
        AMrates(ir), ...
        sum(r1_S1_FRnc)/sum(iRS)*100, sum(r2_S0_FRnc)/sum(iRS)*100,...
        sum(r3_S1_FR0)/sum(iRS)*100,  sum(r4_S0_FR0)/sum(iRS)*100,...
        sum(r5_S1_FR1)/sum(iRS)*100,  sum(r6_S0_FR1)/sum(iRS)*100);
end


% NS cells
fprintf(fileID,'\n----- NS cells -----\n');
for ir = 1:5
    
    % Sync
    r1_S1_FRnc = FR_sig05(iNS,ir)==0  & VS_sig(iNS,ir)==1;
    r3_S1_FR0  = FR_sig05(iNS,ir)==-1 & VS_sig(iNS,ir)==1;
    r5_S1_FR1  = FR_sig05(iNS,ir)==1  & VS_sig(iNS,ir)==1;
    
    % Non-sync
    r2_S0_FRnc = FR_sig05(iNS,ir)==0  & VS_sig(iNS,ir)==0;
    r4_S0_FR0  = FR_sig05(iNS,ir)==-1 & VS_sig(iNS,ir)==0;
    r6_S0_FR1  = FR_sig05(iNS,ir)==1  & VS_sig(iNS,ir)==0;
    
    nsync  = sum(VS_sig(iNS,ir)==1);
    nNsync = sum(VS_sig(iNS,ir)==0);
    vertmid = round(100*sum(VS_sig(iNS,ir)==1)/sum(iNS),1);
    
    subplot(2,5,ir+5)
    xlim([0 100])
    ylim([0 100])
    hold on
    axis square
    set(gca,'Color','none')
    box off
    
    % Sync
    r1_x = [0 vertmid vertmid 0];
    r1_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1) 100 100];
    fill(r1_x,r1_y,0.6.*[1 1 1],'EdgeColor','none')
    
    r3_x = [0 vertmid vertmid 0];
    r3_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)...
        100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1)];
    fill(r3_x,r3_y,[0.9 0.2 0.2],'EdgeColor','none')
    
    r5_x = [0 vertmid vertmid 0];
    r5_y = [0 0 ...
        100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)];
    fill(r5_x,r5_y,[0.2 0.2 0.9],'EdgeColor','none')
    
    % Non-sync
    r2_x = [vertmid 100 100 vertmid];
    r2_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100 100];
    fill(r2_x,r2_y,0*[1 1 1],'EdgeColor','none')
    
    r4_x = [vertmid 100 100 vertmid];
    r4_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1)];
    fill(r4_x,r4_y,[0.6 0 0],'EdgeColor','none')
    
    r6_x = [vertmid 100 100 vertmid];
    r6_y = [0 0 ...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)];
    fill(r6_x,r6_y,[0 0 0.6],'EdgeColor','none')
    
    
    % Print results
    fprintf(fileID,'%i Hz: \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n',...
        AMrates(ir), ...
        sum(r1_S1_FRnc)/sum(iNS)*100, sum(r2_S0_FRnc)/sum(iNS)*100,...
        sum(r3_S1_FR0)/sum(iNS)*100,  sum(r4_S0_FR0)/sum(iNS)*100,...
        sum(r5_S1_FR1)/sum(iNS)*100,  sum(r6_S0_FR1)/sum(iNS)*100);
end

fclose(fileID);

print_eps_kp(gcf,fullfile(savedir,'Overlap_SyncFR_RSNS'))


%%   PIE -- Overlap of FR and Sync responses

figure;
set(gcf,'Position',fullscreen)

% fileID = fopen(fullfile(savedir,'Overlap_SyncFR_RSNS.txt'),'w');

% RS cells
% fprintf(fileID,'\n----- RS cells -----\n');
for ir = 1:5
    
    % Sync
    r1_S1_FRnc = FR_sig05(iRS,ir)==0  & VS_sig(iRS,ir)==1;
    r3_S1_FR0  = FR_sig05(iRS,ir)==-1 & VS_sig(iRS,ir)==1;
    r5_S1_FR1  = FR_sig05(iRS,ir)==1  & VS_sig(iRS,ir)==1;
    
    % Non-sync
    r2_S0_FRnc = FR_sig05(iRS,ir)==0  & VS_sig(iRS,ir)==0;
    r4_S0_FR0  = FR_sig05(iRS,ir)==-1 & VS_sig(iRS,ir)==0;
    r6_S0_FR1  = FR_sig05(iRS,ir)==1  & VS_sig(iRS,ir)==0;
    
    nsync  = sum(VS_sig(iRS,ir)==1)
    nNsync = sum(VS_sig(iRS,ir)==0)
    vertmid = round(100*sum(VS_sig(iRS,ir)==1)/sum(iRS),1);
    
    subplot(4,5,ir)
    set(gca,'Color','none')
    box off
    
    pie([sum(r2_S0_FRnc) sum(r4_S0_FR0) sum(r6_S0_FR1)]./nNsync);
%     [sum(r1_S1_FRnc) sum(r3_S1_FR0) sum(r5_S1_FR1) sum(r2_S0_FRnc) sum(r4_S0_FR0) sum(r6_S0_FR1)]./sum(iRS)
    colormap([0*[1 1 1];...
             [0.6 0 0];...
             [0 0 0.6]])
          
    subplot(4,5,ir+5)
    set(gca,'Color','none')
    box off
    pie([sum(r1_S1_FRnc) sum(r3_S1_FR0) sum(r5_S1_FR1)]./nsync);
    colormap([0*[1 1 1];...
             [0.6 0 0];...
             [0 0 0.6]])
end

print_eps_kp(gcf,fullfile(savedir,'Pie_SyncFR_RSNS'))



% NS cells
% fprintf(fileID,'\n----- NS cells -----\n');
for ir = 1:5
    
    % Sync
    r1_S1_FRnc = FR_sig05(iNS,ir)==0  & VS_sig(iNS,ir)==1;
    r3_S1_FR0  = FR_sig05(iNS,ir)==-1 & VS_sig(iNS,ir)==1;
    r5_S1_FR1  = FR_sig05(iNS,ir)==1  & VS_sig(iNS,ir)==1;
    
    % Non-sync
    r2_S0_FRnc = FR_sig05(iNS,ir)==0  & VS_sig(iNS,ir)==0;
    r4_S0_FR0  = FR_sig05(iNS,ir)==-1 & VS_sig(iNS,ir)==0;
    r6_S0_FR1  = FR_sig05(iNS,ir)==1  & VS_sig(iNS,ir)==0;
    
    nsync  = sum(VS_sig(iNS,ir)==1);
    nNsync = sum(VS_sig(iNS,ir)==0);
    vertmid = round(100*sum(VS_sig(iNS,ir)==1)/sum(iNS),1);
    
    subplot(2,5,ir+5)
    xlim([0 100])
    ylim([0 100])
    hold on
    axis square
    set(gca,'Color','none')
    box off
    
    % Sync
    r1_x = [0 vertmid vertmid 0];
    r1_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1) 100 100];
    fill(r1_x,r1_y,0.6.*[1 1 1],'EdgeColor','none')
    
    r3_x = [0 vertmid vertmid 0];
    r3_y = [100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)...
        100-round(100*sum(r1_S1_FRnc)/nsync,1) 100-round(100*sum(r1_S1_FRnc)/nsync,1)];
    fill(r3_x,r3_y,[0.9 0.2 0.2],'EdgeColor','none')
    
    r5_x = [0 vertmid vertmid 0];
    r5_y = [0 0 ...
        100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)  100-round(100*sum(r1_S1_FRnc)/nsync,1)-round(100*sum(r3_S1_FR0)/nsync,1)];
    fill(r5_x,r5_y,[0.2 0.2 0.9],'EdgeColor','none')
    
    % Non-sync
    r2_x = [vertmid 100 100 vertmid];
    r2_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100 100];
    fill(r2_x,r2_y,0*[1 1 1],'EdgeColor','none')
    
    r4_x = [vertmid 100 100 vertmid];
    r4_y = [100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1) 100-round(100*sum(r2_S0_FRnc)/nNsync,1)];
    fill(r4_x,r4_y,[0.6 0 0],'EdgeColor','none')
    
    r6_x = [vertmid 100 100 vertmid];
    r6_y = [0 0 ...
        100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)  100-round(100*sum(r2_S0_FRnc)/nNsync,1)-round(100*sum(r4_S0_FR0)/nNsync,1)];
    fill(r6_x,r6_y,[0 0 0.6],'EdgeColor','none')
    
    
    % Print results
%     fprintf(fileID,'%i Hz: \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n    %0.1f%%    %0.1f%% \n',...
%         AMrates(ir), ...
%         sum(r1_S1_FRnc)/sum(iNS)*100, sum(r2_S0_FRnc)/sum(iNS)*100,...
%         sum(r3_S1_FR0)/sum(iNS)*100,  sum(r4_S0_FR0)/sum(iNS)*100,...
%         sum(r5_S1_FR1)/sum(iNS)*100,  sum(r6_S0_FR1)/sum(iNS)*100);
end

% fclose(fileID);

print_eps_kp(gcf,fullfile(savedir,'Pie_SyncFR_RSNS'))



%%
keyboard
% 
% % FF: avg per stimulus
% FFplot = nan(1,5);
% for ist = 1:length(FFplot)
%     FFplot(ist) = median(FFmeans(~isnan(FFmeans(:,ist+1)),ist+1));
% end
% maxFF=max(max(FFmeans(:,2:end)));
% 
% subplot(3,4,10)
% plotSpread(FFmeans(:,2:end),'distributionColors','k')
% plot(1:5,FFplot,'.-','MarkerSize',50,'LineWidth',4,'Color',[0.9 0.8 0])
% set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
% xlim([-1 6])
% ylim([0 10])
% xlabel('AM rate (Hz)')
% ylabel('Mean FF')
% axis square
% box on
% 
% 
% % FR range
% % [~,iu]=sort(FRmaxs(:,1)-FRmins(:,1));
% [~,iu]=sort([UnitData.BaseFR]);
% subplot(3,4,[1 5 9]);
% for ir=1:5
%     plot([FRmins(iu,ir) FRmaxs(iu,ir)]',ir+0.8.*[1:iUn; 1:iUn]./iUn-0.4,'-k')
%     hold on
% end
% % plot(iu,FFplot,'.-','MarkerSize',50,'LineWidth',4,'Color',[0.9 0.8 0])
% set(gca,'ytick',1:5,'yticklabel',AMrates,'xtick',0:20:80,'xticklabel',0:20:80,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
% ylim([0.5 5.5])
% xlim([0 80])
% ylabel('AM rate (Hz)')
% xlabel('FR min to max')
% box on
% axis fill



% Distribution of baseline-subtracted FRs per stimulus
% % hold off; cla
% % window=gausswin(5);
% % window=window-min(window);
% % window=window/sum(window);
% % 
% % for ist = 1:size(dFRs,2)
% % %     ih(ist) = histogram(dFRs(:,ist),-50:0.1:100,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
% %     [values,edges] = histcounts(dFRs(:,ist),-31:0.5:31);
% %     sig=conv(values,window,'same');
% %     ih(ist) = plot(edges(2:end)-mode(diff(edges))/2,sig./sum(sig),'Color',colors(ist+1,:),'LineWidth',2);
% %     hold on
% % end
% % axis square
% % set(gca,'Color','none')
% % xlim(30*[-1 1])
% % ylim([0 0.1])
% % xlabel('Sound driven FR')
% % ylabel('Probability')
% % legend(ih,Info.stim_ID_key(2:end),'FontSize',9,'Location','northeast')
% % legend('boxoff')
% % grid on


p_F  = nan(1,5);
p_A  = nan(1,5);
p_NC = nan(1,5);
p_ns = nan(1,5);
for irate = 1:5
    p_F(irate)  = sum(strcmp(PdcRespType(:,irate),'F')) / (size(PdcRespType,1)-sum(strcmp(PdcRespType(:,irate),'n')));
    p_A(irate)  = sum(strcmp(PdcRespType(:,irate),'A')) / (size(PdcRespType,1)-sum(strcmp(PdcRespType(:,irate),'n')));
    p_NC(irate) = sum(strcmp(PdcRespType(:,irate),'S')) / (size(PdcRespType,1)-sum(strcmp(PdcRespType(:,irate),'n')));
    p_ns(irate) = sum(strcmp(PdcRespType(:,irate),'n')) / size(PdcRespType,1);
end

subplot(3,4,1);
cla
plot(round(p_F*100,1),'g-','LineWidth',3)
hold on
plot(round(p_A*100,1),'r-','LineWidth',3)
% plot(p_NC,'k-','LineWidth',3)
% plot(p_ns,'-','Color',0.5*[1 1 1],'LineWidth',3)

axis square
ylim([0 25])
xlim([0 6])
set(gca,'xtick',1:5,'xticklabel',AMrates,'Color','none')
ylabel('Percent of cells')
legend({'F' 'A'},'Location','northwest')




% Save figures

savedir = fullfile(fn.figs,'PopulationTuning','Apr2020');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,'TuningSummary_RSNS'))


end





