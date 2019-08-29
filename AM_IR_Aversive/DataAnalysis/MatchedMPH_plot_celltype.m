function MatchedMPH_plot_celltype(USE_MEASURE)
%  MPH_matchedComparison
%
%   Loads data from MPH_matchedComparison. 
%
%   Plots grouped by cell type/spike waveform shape.
%
% KP, 2019-08
% 

global fn  


%% Load data 

fn = set_paths_directories('','',1);

% Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

if nargin<1
    USE_MEASURE = 'FR';
end

switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 0.01;
                axmax = 15*ceil(max([UnitData.BaseFR])/10);
            case 'z'
                axmin = -1;
                axmax = 1.5;
        end
    case {'TrV' 'FF'}
        units = ' ';
        axmin = 0.01;
        axmax = 4;
end


% Matched MP data

savedir = fullfile(fn.figs,'MPHmatched');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = ['matchedMPH_' USE_MEASURE];

% Load already saved data
load(fullfile(savedir,savename))

% Change savename for figure
savename = ['matchedMPH_' USE_MEASURE '_CellType'];


%%
% Filter data by AM rate

% irate = 5;
% 
% PD_PdcIR     = PD_PdcIR(MPrate==irate,:);
% PD_PdcPdc    = PD_PdcPdc(MPrate==irate,:);
% % PD_IRIR      = PD_IRIR(MPrate==irate,:);
% PdcRespType  = PdcRespType(MPrate==irate);
% PdcStartTime = PdcStartTime(MPrate==irate);
% 
% figure; 
% plot(diff(PD_PdcPdc),diff(PD_PdcIR),'ko')


%% SET UP FIGURES
%::::::::::::::::::

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

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

PlotMarkerSize = 15;

histbinedges = linspace(-sqrt(60^2 + 60^2)/2,sqrt(60^2 + 60^2)/2,40);


%::::::::::::::::::
% CREATE PLOTS 
%::::::::::::::::::
hf = figure;
set(hf,'Position',fullscreen,'NextPlot','add')

%---------------
% Broad spikes
%---------------
% Pdc-Irr
hs(1)=subplot(6,4,[1 5 9]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-IR')

% Pdc-Pdc
hs(2)=subplot(6,4,[2 6 10]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-Pdc')

% DIFFERENCE HISTOGRAM
hs(7)=subplot(6,4,[4 8]);
hold on
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
set(gca,'Color','none','xtick',[],'ytick',[])
title('Broad spikes')


%---------------
% Narrow spikes
%---------------
% Pdc-IR 
hs(4)=subplot(6,4,[13 17 21]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-IR')

% Pdc-Pdc
hs(5)=subplot(6,4,[14 18 22]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-Pdc')

% DIFFERENCE HISTOGRAM
hs(8)=subplot(6,4,[12 16]);
hold on
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
set(gca,'Color','none','xtick',[],'ytick',[])
title('Narrow spikes')


% suptitle([num2str(AMrates(irate)) ' Hz'])

%--------------- blank for now
% Irr-Irr  
% hs(3)=subplot(6,4,[3 7 11]);
% hold on
% plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
% set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
% axis square
% title('IR-IR')
% 
% % IR-IR density plot   OR   starttimeXdiff 
% hs(6)=subplot(6,4,[15 19 23]);
% set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
% axis square
% title('IR-IR')



%% PLOT DATA  
%:::::::::::::::::::

iB = strcmp(CellType,'B');
iN = strcmp(CellType,'N');

%---------------
% Broad spikes
%---------------

%--Pdc-Irr
subplot(hs(1)); hold on
plot(PD_PdcIR(iB,1),  PD_PdcIR(iB,2), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.99*[0.9 0.2 0.2],'MarkerEdgeColor','none');

%--Pdc-Pdc
subplot(hs(2)); hold on
plot(PD_PdcPdc(iB,1), PD_PdcPdc(iB,2),'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.5*[0.9 0.2 0.2],'MarkerEdgeColor','none');

%--Diff hist
subplot(hs(7)); hold on
histogram(PD_PdcIR(iB,1)-PD_PdcIR(iB,2),  histbinedges,'Normalization','pdf','FaceColor',0.99*[0.9 0.2 0.2],'EdgeColor','none','FaceAlpha',1);
histogram(PD_PdcPdc(iB,1)-PD_PdcPdc(iB,2),histbinedges,'Normalization','pdf','FaceColor','none','EdgeColor',0.2*[0.9 0.2 0.2],'FaceAlpha',1,'DisplayStyle','stairs','LineWidth',2);
plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
ylim([0 0.3])


%---------------
% Narrow spikes
%---------------

%--Pdc-Irr
subplot(hs(4)); hold on
plot(PD_PdcIR(iN,1),  PD_PdcIR(iN,2), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.99*[0.2 0.2 0.9],'MarkerEdgeColor','none');

%--Pdc-Pdc
subplot(hs(5)); hold on
plot(PD_PdcPdc(iN,1), PD_PdcPdc(iN,2),'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.5*[0.2 0.2 0.9],'MarkerEdgeColor','none');

%--Diff hist
subplot(hs(8)); hold on
histogram(PD_PdcIR(iN,1)-PD_PdcIR(iN,2),  histbinedges,'Normalization','pdf','FaceColor',0.99*[0.2 0.2 0.9],'EdgeColor','none','FaceAlpha',1);
histogram(PD_PdcPdc(iN,1)-PD_PdcPdc(iN,2),histbinedges,'Normalization','pdf','FaceColor','none','EdgeColor',0.2*[0.2 0.2 0.9],'FaceAlpha',1,'DisplayStyle','stairs','LineWidth',2);
plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
ylim([0 0.3])


%% STATISTICS  
%:::::::::::::::::


%---------------
%   Pdc-Irr
%---------------

%--Broad spikes
[rc,pc] = corrcoef(PD_PdcIR(iB,1), PD_PdcIR(iB,2));
m   = PD_PdcIR(iB,1)\PD_PdcIR(iB,2);
pwsr    = signrank(PD_PdcIR(iB,1), PD_PdcIR(iB,2));
pwrs    = ranksum(PD_PdcIR(iB,1),  PD_PdcIR(iB,2));

subplot(hs(1)); hold on
stattext = sprintf('WSR p=%0.2f\nr=%0.2f, m=%0.2f',pwsr,rc(1,2),m);
text(axmax/4*3,axmax/10,stattext)

%--Narrow spikes
[rc,pc] = corrcoef(PD_PdcIR(iN,1), PD_PdcIR(iN,2));
m   = PD_PdcIR(iN,1)\PD_PdcIR(iN,2);
pwsr    = signrank(PD_PdcIR(iN,1), PD_PdcIR(iN,2));
pwrs    = ranksum(PD_PdcIR(iN,1),  PD_PdcIR(iN,2));

subplot(hs(4)); hold on
stattext = sprintf('WSR p=%0.2f\nr=%0.2f, m=%0.2f',pwsr,rc(1,2),m);
text(axmax/4*3,axmax/10,stattext)


%---------------
%   Pdc-Pdc
%---------------

%--Broad spikes
[rc,pc] = corrcoef(PD_PdcPdc(iB,1), PD_PdcPdc(iB,2));
m       = PD_PdcPdc(iB,1)\PD_PdcPdc(iB,2);
pwsr    = signrank(PD_PdcPdc(iB,1), PD_PdcPdc(iB,2));
pwrs    = ranksum(PD_PdcPdc(iB,1),  PD_PdcPdc(iB,2));

subplot(hs(2)); hold on
stattext = sprintf('WSR p=%0.2f\nr=%0.2f, m=%0.2f',pwsr,rc(1,2),m);
text(axmax/4*3,axmax/10,stattext)

%--Narrow spikes
[rc,pc] = corrcoef(PD_PdcPdc(iN,1), PD_PdcPdc(iN,2));
m   = PD_PdcPdc(iN,1)\PD_PdcPdc(iN,2);
pwsr    = signrank(PD_PdcPdc(iN,1), PD_PdcPdc(iN,2));
pwrs    = ranksum(PD_PdcPdc(iN,1),  PD_PdcPdc(iN,2));

subplot(hs(5)); hold on
stattext = sprintf('WSR p=%0.2f\nr=%0.2f, m=%0.2f',pwsr,rc(1,2),m);
text(axmax/4*3,axmax/10,stattext)



%~~~~~~~~~~~~~~
% DIFFERENCES
%~~~~~~~~~~~~~~

%---------------
% Broad spikes
%---------------

PdcDiffs  = PD_PdcPdc(iB,1)-PD_PdcPdc(iB,2);
PdIRDiffs = PD_PdcIR(iB,1)-PD_PdcIR(iB,2);

% Normality test
clear h_l
h_l(1) = lillietest(PdcDiffs);  % h=1 means distr NOT normal
h_l(2) = lillietest(PdIRDiffs);
if any(h_l==1)% non-parametric
    p_rs       = ranksum(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','BrownForsythe','display','off');
    stattext   = sprintf('ranksum p=%0.2f \nBrownForsythe p=%0.2e',p_rs,p_var);
else          % parametric
    [~,p_t]    = ttest2(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','Bartlett','display','off');
    [~,p_var]  = vartest2(PdcDiffs, PdIRDiffs);
    stattext   = sprintf('ttest p=%0.2f \nF test p=%0.2e',p_t,p_var);
    var_PdcPdc = var(PdcDiffs);
    var_PdcIR  = var(PdIRDiffs);
end

subplot(hs(7)); hold on
text(-40,0.28,stattext)


%---------------
% Narrow spikes
%---------------

PdcDiffs  = PD_PdcPdc(iN,1)-PD_PdcPdc(iN,2);
PdIRDiffs = PD_PdcIR(iN,1)-PD_PdcIR(iN,2);

% Normality test
clear h_l
h_l(1) = lillietest(PdcDiffs);  % h=1 means distr NOT normal
h_l(2) = lillietest(PdIRDiffs);
if any(h_l==1)% non-parametric
    p_rs       = ranksum(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','BrownForsythe','display','off');
    stattext   = sprintf('ranksum p=%0.2f \nBrownForsythe p=%0.2e',p_rs,p_var);
else          % parametric
    [~,p_t]    = ttest2(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','Bartlett','display','off');
    [~,p_var]  = vartest2(PdcDiffs, PdIRDiffs);
    stattext   = sprintf('ttest p=%0.2f \nF test p=%0.2e',p_t,p_var);
    var_PdcPdc = var(PdcDiffs);
    var_PdcIR  = var(PdIRDiffs);
end

subplot(hs(8)); hold on
text(-40,0.28,stattext)


%% SAVE FIGURES

print_eps_kp(hf,fullfile(savedir,savename));


end %function




