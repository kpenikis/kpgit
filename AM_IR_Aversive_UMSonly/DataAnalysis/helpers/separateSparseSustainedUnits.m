function separateSparseSustainedUnits(Resp)
%
%  pp_plot_rasters(subject, session, [channel, clu] )
%    Plots a raster and psth for each stimulus, separating trials by
%    the preceding stimulus. 
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2018-05
%



%% Load Resp data table 

fn = set_paths_directories('','',1);

if nargin<1
    q = load(fullfile(fn.processed,'RespStruct_allSU'));
    Resp = q.Resp;
end

%% Prepare figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
tallrect = [1 scrsz(4) scrsz(3)/4 scrsz(4)];

hf = figure;
set(hf,'Position',tallrect)


%% Inspect the distribution of FRs

minTrs   =  10;
histbins = logspace(-1,1.61,20);


yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {Resp.FR_raw_tr},'UniformOutput',false);
yyy = vertcat(yyy{:});


%% Plot the distributions of FR repsonses across units for some stimuli


% All stimuli
AllStim_FRs = mean(yyy,2,'omitnan');
[~,isrt] = sort(AllStim_FRs);
AllStim_FRs = AllStim_FRs(isrt);

isp(1)=subplot(5,1,1);
histogram(AllStim_FRs,histbins,'FaceColor',0.3*[1 1 1])
% set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% xlim([0 max(histbins)])
title('Avg of ALL stim responses')


% Pdc stimuli
PdcStim_FRs = mean(yyy(:,2:6),2,'omitnan');
[~,isrt] = sort(PdcStim_FRs);
PdcStim_FRs = PdcStim_FRs(isrt);

hold on
isp(2)=subplot(5,1,2);
histogram(PdcStim_FRs,histbins,'FaceColor',0.3*[1 1 1])
% set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% xlim([0 max(histbins)])
title('Avg of Periodic stim responses')


% Warn
Warn_FR = yyy(:,1);
[~,isrt] = sort(Warn_FR);
Warn_FR = Warn_FR(isrt);

isp(3)=subplot(5,1,3);
histogram(Warn_FR,histbins,'FaceColor',0.3*[1 1 1])
% set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% xlim([0 max(histbins)])
% ylim([0 ymaxval])
title('Unmodulated noise responses')


% IR A
AC_FR = yyy(:,7);
[~,isrt] = sort(AC_FR);
AC_FR = AC_FR(isrt);

hold on
isp(4)=subplot(5,1,4);
histogram(AC_FR,histbins,'FaceColor',0.3*[1 1 1])
title('Irregular seq A responses')


% IR B
DB_FR = yyy(:,8);
[~,isrt] = sort(DB_FR);
DB_FR = DB_FR(isrt);

hold on
isp(5)=subplot(5,1,5);
histogram(DB_FR,histbins,'FaceColor',0.3*[1 1 1])
% set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
xlabel('Mean FR (spikes/sec)')
ylabel('Count')
title('Irregular seq B responses')

% Finish plots
ymaxval = 18;

linkaxes(isp,'xy')
set(isp,'xscale','log','xtick',[0.1 1 5 10 20 40],'ytick',0:5:ymaxval)
xlim([0 max(histbins)])
ylim([0 ymaxval])


% Save fig
savedir = fullfile(fn.processed,'Rasters_SortFR');
% print_eps_kp(gcf,fullfile(savedir,'Population_FR_Distributions'))
% print_svg_kp(gcf,fullfile(savedir,'Population_FR_Distributions'))


%% Pick a boundary to separate units

% Warn
Warn_FR = yyy(:,1);
BOUNDARY = 3; %Hz

SparseUnits = Resp(Warn_FR<BOUNDARY);
SustUnits   = Resp(Warn_FR>BOUNDARY);


% Periodic
PdcStim_FRs = mean(yyy(:,2:6),2,'omitnan');
BOUNDARY = 3; %Hz

SparseUnits = Resp(PdcStim_FRs<BOUNDARY);
SustUnits   = Resp(PdcStim_FRs>BOUNDARY);


% Average of IR sequences    ***
AC_FR = yyy(:,7);
DB_FR = yyy(:,8);
BOUNDARY = 5; %Hz

SparseUnits = Resp(nanmean([DB_FR AC_FR],2)<BOUNDARY);
SustUnits   = Resp(nanmean([DB_FR AC_FR],2)>BOUNDARY);


% CHECK OUTPUT AGAINST RASTERS


% After making decision of where and how to draw boundary, add field to
% Resp struct

[Resp(:).Sparse] = deal(0);
[Resp(nanmean([DB_FR AC_FR],2)<BOUNDARY).Sparse] = deal(1);



%% Mean phase distributions

AMrates = [2 4 8 16 32];

MeanPhase_sparse = vertcat(SparseUnits.Phase_deg);
MeanPhase_highFR = vertcat(SustUnits.Phase_deg);

hfmp = figure;
set(hfmp,'Position',tallrect)
for ir = 1:5
    subplot(5,1,ir); hold on
    ih(1,ir)=histogram( mod(MeanPhase_highFR(:,ir),360), 0:36:360,'FaceColor','k');
    ih(2,ir)=histogram( mod(MeanPhase_sparse(:,ir),360), 0:36:360,'FaceColor','b');
    xlim([0 360])
    hold off
    title([num2str(AMrates(ir)) ' Hz'])
end
xlabel('Mean Phase (deg)')
ylabel('Number of units')

subplot(5,1,1); hold on
legend(ih(:,1),{'High FR' 'Sparse'})
suptitle(sprintf('All units, all stimuli\nraw mean phase response'))

print_eps_kp(hfmp,fullfile(savedir,'MeanPhase_Dist_raw'))
print_svg_kp(hfmp,fullfile(savedir,'MeanPhase_Dist_raw'))



hfmpa = figure;
set(hfmpa,'Position',tallrect)
for ir = 1:5
    subplot(5,1,ir); hold on
    ih(1,ir)=histogram( mod(MeanPhase_highFR(:,ir) - vertcat(SustUnits.IntTime)./(1000./AMrates).*360 ,360), 0:36:360,'FaceColor','k');
    ih(2,ir)=histogram( mod(MeanPhase_sparse(:,ir) - vertcat(SparseUnits.IntTime)./(1000./AMrates).*360 ,360), 0:36:360,'FaceColor','b');
    xlim([0 360])
    hold off
    title([num2str(AMrates(ir)) ' Hz'])
end
xlabel('Mean Phase (deg)')
ylabel('Number of units')

subplot(5,1,1); hold on
legend(ih(:,1),{'High FR' 'Sparse'})
suptitle(sprintf('All units, all stimuli\nmean phase adjusted acc. to each unit''s Integration Time'))

% print_eps_kp(hfmpa,fullfile(savedir,'MeanPhase_Dist_adjIntTime'))
% print_svg_kp(hfmpa,fullfile(savedir,'MeanPhase_Dist_adjIntTime'))



%% Mean phase distributions, only synchronized rates

MeanPhase_sparse = nan(size(MeanPhase_sparse));
[~,SparseUnits] = identifyResponsiveUnits(SparseUnits);
for iun = 1:numel(SparseUnits)
    idx = SparseUnits(iun).iSync;                           % to adjust phase based on unit's integration time
    MeanPhase_sparse(iun,idx) = [SparseUnits(iun).Phase_deg(idx)] - SparseUnits(iun).IntTime./(1000./AMrates(idx)).*360;
end

[~,SustUnits] = identifyResponsiveUnits(SustUnits);
MeanPhase_highFR = nan(size(MeanPhase_highFR));
for iun = 1:numel(SustUnits)
    idx = SustUnits(iun).iSync;                             % to adjust phase based on unit's integration time
    MeanPhase_highFR(iun,idx) = [SustUnits(iun).Phase_deg(idx)] - SustUnits(iun).IntTime./(1000./AMrates(idx)).*360;
end


hfmp = figure;
set(hfmp,'Position',tallrect)

for ir = 1:5
    
    subplot(5,1,ir); hold on
    ih(1,ir)=histogram( mod(MeanPhase_highFR(:,ir),360), 0:36:360,'FaceColor','k');
    ih(2,ir)=histogram( mod(MeanPhase_sparse(:,ir),360), 0:36:360,'FaceColor','b');
    xlim([0 360])
    hold off
    title([num2str(AMrates(ir)) ' Hz'])
    
end

xlabel('Mean Phase (deg)')
ylabel('Number of units')
% suptitle(sprintf('All units, Sync stimuli only\nraw mean phase response'))
suptitle(sprintf('All units, Sync stimuli only\nmean phase adjusted acc. to each unit''s Integration Time'))

legend(ih(:,1),{'High FR' 'Sparse'})

% print_eps_kp(gcf,fullfile(savedir,'MeanPhase_Dist_Sync_adjIntTime'))
% print_svg_kp(gcf,fullfile(savedir,'MeanPhase_Dist_Sync_adjIntTime'))
% 
% print_eps_kp(gcf,fullfile(savedir,'MeanPhase_Dist_Sync_raw'))
% print_svg_kp(gcf,fullfile(savedir,'MeanPhase_Dist_Sync_raw'))



%% Save Resp struct again, with new field

filename = 'RespStruct_allSU';
save(fullfile(fn.processed,filename),'Resp','-v7.3');





end
