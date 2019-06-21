function plotUnitDistributions(UnitData,UnitInfo)
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

close all
global fn


%% Load Unit data files

fn = set_paths_directories('','',1);

if nargin<1
    q = load(fullfile(fn.processed,'Units'));
    UnitData = q.UnitData;
    UnitInfo = q.UnitInfo;
    clear q
end

[~,UnitData] = identifyResponsiveUnits(UnitData);


%% Prepare figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
tallrect = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
widerect = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

hf = figure;
set(hf,'Position',tallrect)


% And set some other params
histbins = logspace(-1,1.61,20);


%% Prepare RespTypes

% % % % % % idxSparse    = find(strcmp(UnitInfo.RespType,'sparse'));
% % % % % % idxSustained = find(strcmp(UnitInfo.RespType,'sustained'));
% % % % % % idxGap       = find(strcmp(UnitInfo.RespType,'gap'));
% % % % % % idxNone      = find(strcmp(UnitInfo.RespType,'na'));
% % % % % 
% % % % % idxSparse    = find(UnitInfo.OWR_2hz>0.2);
% % % % % % idxSustained = find(strcmp(UnitInfo.RespType,'sustained'));
% % % % % idxGap       = find(UnitInfo.OWR_2hz<-0.2);
% % % % % idxNone      = find(UnitInfo.OWR_2hz>-0.2 & UnitInfo.OWR_2hz<0.2);
% % % % % idxSustained      = find(UnitInfo.OWR_2hz>-0.2 & UnitInfo.OWR_2hz<0.2);
% % % % % 
% % % % % 
% % % % % 
% % % % % colSparse    = [0.157 0.44  1]; %[0.208 0.678 0.275];
% % % % % % colSustained = [1     0.44  0.157]; %[0.200 0.282 0.800];
% % % % % colGap       = [0.139 0.159 0.190];
% % % % % colSustained = 0.7 * [1 1 1];
% % % % % 
% % % % % 
% % % % % alphaSparse     = 0.3;
% % % % % alphaSustained  = 0.7;
% % % % % alphaGap        = 0.8;
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% FR
% % % % % 
% % % % % AllStimResponses = cellfun(@(xxx) mean(xxx,1,'omitnan'), {UnitData.FR_raw_tr},'UniformOutput',false);
% % % % % AllStimResponses = vertcat(AllStimResponses{:});
% % % % % 
% % % % % AllBaselineFRs = [UnitData.BaseFR]';
% % % % % 
% % % % % 
% % % % % %% Plot the distributions of FR repsonses across units for some stimuli
% % % % % 
% % % % % % All stimuli
% % % % % AllStim_FRs = mean(AllStimResponses,2,'omitnan');
% % % % % 
% % % % % isp(1)=subplot(5,1,1); hold on
% % % % % % histogram(AllStim_FRs(idxNone),histbins,'FaceColor',colNone,'FaceAlpha',0.4)
% % % % % histogram(AllStim_FRs(idxGap),histbins,'FaceColor',colGap,'FaceAlpha',alphaGap)
% % % % % histogram(AllStim_FRs(idxSustained),histbins,'FaceColor',colSustained,'FaceAlpha',alphaSustained)
% % % % % histogram(AllStim_FRs(idxSparse),histbins,'FaceColor',colSparse,'FaceAlpha',alphaSparse)
% % % % % set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% % % % % xlim([0 max(histbins)])
% % % % % title('Avg sound response')
% % % % % 
% % % % % 
% % % % % % Pdc stimuli
% % % % % PdcStim_FRs = mean(AllStimResponses(:,2:6),2,'omitnan');
% % % % % 
% % % % % isp(2)=subplot(5,1,2); hold on
% % % % % % histogram(PdcStim_FRs(idxNone),histbins,'FaceColor',colNone,'FaceAlpha',0.4)
% % % % % histogram(PdcStim_FRs(idxGap),histbins,'FaceColor',colGap,'FaceAlpha',alphaGap)
% % % % % histogram(PdcStim_FRs(idxSustained),histbins,'FaceColor',colSustained,'FaceAlpha',alphaSustained)
% % % % % histogram(PdcStim_FRs(idxSparse),histbins,'FaceColor',colSparse,'FaceAlpha',alphaSparse)
% % % % % set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% % % % % xlim([0 max(histbins)])
% % % % % title('Avg of Periodic stim responses')
% % % % % 
% % % % % 
% % % % % % Avg IR
% % % % % IR_FR = nanmean(AllStimResponses(:,7:8),2);
% % % % % 
% % % % % isp(3)=subplot(5,1,3); hold on
% % % % % % histogram(IR_FR(idxNone),histbins,'FaceColor',colNone)
% % % % % histogram(IR_FR(idxGap),histbins,'FaceColor',colGap,'FaceAlpha',alphaGap)
% % % % % histogram(IR_FR(idxSustained),histbins,'FaceColor',colSustained,'FaceAlpha',alphaSustained)
% % % % % histogram(IR_FR(idxSparse),histbins,'FaceColor',colSparse,'FaceAlpha',alphaSparse)
% % % % % set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% % % % % xlim([0 max(histbins)])
% % % % % title('Avg Irregular seq responses')
% % % % % 
% % % % % 
% % % % % % Warn
% % % % % Warn_FR = AllStimResponses(:,1);
% % % % % 
% % % % % isp(4)=subplot(5,1,4); hold on
% % % % % % histogram(Warn_FR(idxNone),histbins,'FaceColor',colNone,'FaceAlpha',0.4)
% % % % % histogram(Warn_FR(idxGap),histbins,'FaceColor',colGap,'FaceAlpha',alphaGap)
% % % % % histogram(Warn_FR(idxSustained),histbins,'FaceColor',colSustained,'FaceAlpha',alphaSustained)
% % % % % histogram(Warn_FR(idxSparse),histbins,'FaceColor',colSparse,'FaceAlpha',alphaSparse)
% % % % % set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% % % % % xlim([0 max(histbins)])
% % % % % title('Unmodulated noise (Warn) responses')
% % % % % 
% % % % % 
% % % % % % Silence / Baseline FR
% % % % % 
% % % % % isp(5)=subplot(5,1,5); hold on
% % % % % % histogram(AllBaselineFRs(idxNone),histbins,'FaceColor',colNone,'FaceAlpha',0.4)
% % % % % histogram(AllBaselineFRs(idxGap),histbins,'FaceColor',colGap,'FaceAlpha',alphaGap)
% % % % % histogram(AllBaselineFRs(idxSustained),histbins,'FaceColor',colSustained,'FaceAlpha',alphaSustained)
% % % % % histogram(AllBaselineFRs(idxSparse),histbins,'FaceColor',colSparse,'FaceAlpha',alphaSparse)
% % % % % set(gca,'xscale','log','xtick',[0.1 1 5 10 20 40])
% % % % % xlim([0 max(histbins)])
% % % % % xlabel('Mean FR (spikes/sec)')
% % % % % ylabel('Count')
% % % % % title('Baseline FRs (silence)')
% % % % % 
% % % % % 
% % % % % % Finish plots
% % % % % ymaxval = 9;
% % % % % 
% % % % % linkaxes(isp,'xy')
% % % % % set(isp,'xscale','log','xtick',[0.1 1 5 10 20 40],'ytick',0:5:ymaxval)
% % % % % xlim([0 max(histbins)])
% % % % % ylim([0 ymaxval])
% % % % % 
% % % % % 
% % % % % % Save fig
% % % % % savedir = fullfile(fn.processed,'Units','OWRcategorized');
% % % % % if ~exist(savedir,'dir')
% % % % %     mkdir(savedir)
% % % % % end
% % % % % 
% % % % % print_eps_kp(gcf,fullfile(savedir,'FR_Distributions'))
% % % % % print_svg_kp(gcf,fullfile(savedir,'FR_Distributions'))
% % % % % 
% % % % % 
% % % % % % Plot sound FR as function of baseline FR
% % % % % figure;
% % % % % plot([0.1 40],[0.1 40],'--k')
% % % % % hold on
% % % % % plot(AllBaselineFRs(idxGap),AllStim_FRs(idxGap),'o','Color',colGap,'LineWidth',2,'MarkerSize',10)
% % % % % plot(AllBaselineFRs(idxSustained),AllStim_FRs(idxSustained),'o','Color',colSustained,'LineWidth',2,'MarkerSize',10)
% % % % % plot(AllBaselineFRs(idxSparse),AllStim_FRs(idxSparse),'o','Color',colSparse,'LineWidth',2,'MarkerSize',10)
% % % % % axis square
% % % % % set(gca,'xscale','log','yscale','log',...
% % % % %     'xlim',[0.1 40],'ylim',[0.1 40],...
% % % % %     'xtick',[0.1 1 5 10 20 40],'ytick',[0.1 1 5 10 20 40])
% % % % % xlabel('Baseline FR')
% % % % % ylabel('All sound avg FR')
% % % % % 
% % % % % % print_eps_kp(gcf,fullfile(savedir,'BaselineSoundFR'))
% % % % % % print_svg_kp(gcf,fullfile(savedir,'BaselineSoundFR'))
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% Store data in UnitInfo
% % % % % UnitInfo.BaselineFR  = AllBaselineFRs;
% % % % % UnitInfo.Warn_FR     = Warn_FR;
% % % % % UnitInfo.IR_FR       = IR_FR;
% % % % % UnitInfo.PdcStim_FRs = PdcStim_FRs;
% % % % % UnitInfo.AllStim_FRs = AllStim_FRs;
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% Mean phase
% % % % % 
% % % % % AMrates = [2 4 8 16 32];
% % % % % 
% % % % % MeanPhase_Sparse    = vertcat(UnitData(idxSparse).Phase_spk);
% % % % % MeanPhase_Sustained = vertcat(UnitData(idxSustained).Phase_spk);
% % % % % MeanPhase_Gap       = vertcat(UnitData(idxGap).Phase_spk);
% % % % % 
% % % % % hfmp = figure;
% % % % % set(hfmp,'Position',tallrect)
% % % % % for ir = 1:5
% % % % %     isp(ir)=subplot(5,1,ir); hold on
% % % % %     ih(1,ir)=histogram( mod(MeanPhase_Gap(:,ir),360), 0:36:360,'FaceColor',colGap,'FaceAlpha',alphaGap);
% % % % %     ih(2,ir)=histogram( mod(MeanPhase_Sustained(:,ir),360), 0:36:360,'FaceColor',colSustained,'FaceAlpha',alphaSustained);
% % % % %     ih(3,ir)=histogram( mod(MeanPhase_Sparse(:,ir),360), 0:36:360,'FaceColor',colSparse,'FaceAlpha',alphaSparse);
% % % % %     xlim([0 360])
% % % % %     hold off
% % % % %     title([num2str(AMrates(ir)) ' Hz'])
% % % % % end
% % % % % xlabel('Mean Phase (deg)')
% % % % % ylabel('Number of units')
% % % % % 
% % % % % ymaxval = 12;
% % % % % linkaxes(isp,'xy')
% % % % % ylim([0 ymaxval])
% % % % % 
% % % % % % subplot(5,1,1); hold on
% % % % % legend(ih(:,5),{'Gap' 'Sustained' 'Sparse'},'Location','best')
% % % % % suptitle(sprintf('ex units, all stimuli\nmean phase response'))
% % % % % 
% % % % % % print_eps_kp(hfmp,fullfile(savedir,'MeanPhaseDist_allStim'))
% % % % % % print_svg_kp(hfmp,fullfile(savedir,'MeanPhaseDist_allStim'))
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% Spk phase vs. Gap phase
% % % % % 
% % % % % MeanPhase_Sparse    = vertcat(UnitData(idxSparse).Phase_spk);
% % % % % MeanPhase_Sustained = vertcat(UnitData(idxSustained).Phase_spk);
% % % % % MeanPhase_Gap       = vertcat(UnitData(idxGap).Phase_spk);
% % % % % 
% % % % % MeanPhase_gap_Sparse    = vertcat(UnitData(idxSparse).Phase_gap);
% % % % % MeanPhase_gap_Sustained = vertcat(UnitData(idxSustained).Phase_gap);
% % % % % MeanPhase_gap_Gap       = vertcat(UnitData(idxGap).Phase_gap);
% % % % % 
% % % % % 
% % % % % hfsg = figure;
% % % % % set(hfsg,'Position',widerect)
% % % % % for ir = 1:5
% % % % %     subplot(1,5,ir); 
% % % % %     hold on
% % % % %     ih(1,ir)=plot( mod(MeanPhase_Gap(:,ir),360), mod(MeanPhase_gap_Gap(:,ir),360), 'o', 'Color',colGap,'LineWidth',2,'MarkerSize',10);
% % % % %     ih(2,ir)=plot( mod(MeanPhase_Sustained(:,ir),360), mod(MeanPhase_gap_Sustained(:,ir),360), 'o', 'Color',colSustained,'LineWidth',2,'MarkerSize',10);
% % % % %     ih(3,ir)=plot( mod(MeanPhase_Sparse(:,ir),360), mod(MeanPhase_gap_Sparse(:,ir),360), 'o', 'Color',colSparse,'LineWidth',2,'MarkerSize',10);
% % % % %     xlim([0 360])
% % % % %     ylim([0 360])
% % % % %     xlabel('Mean Phase (spikes)')
% % % % %     ylabel('Mean Phase (gaps)')
% % % % %     axis square
% % % % %     hold off
% % % % %     title([num2str(AMrates(ir)) ' Hz'])
% % % % % end
% % % % % 
% % % % % % print_eps_kp(hfsg,fullfile(savedir,'MeanPhase_SpkGap'))
% % % % % % print_svg_kp(hfsg,fullfile(savedir,'MeanPhase_SpkGap'))
% % % % % 
% % % % % 
% % % % % 
% % % % % %% Mean phase distributions, only synchronized rates
% % % % % 
% % % % % MeanPhase_Sparse = nan(length(idxSparse),5);
% % % % % for iun = idxSparse'
% % % % %     idx = UnitData(iun).iSync;              
% % % % %     MeanPhase_Sparse(iun==idxSparse,idx) = [UnitData(iun).Phase_spk(idx)];
% % % % % end
% % % % % 
% % % % % MeanPhase_Sustained = nan(length(idxSustained),5);
% % % % % for iun = idxSustained'
% % % % %     idx = UnitData(iun).iSync;              
% % % % %     MeanPhase_Sustained(iun==idxSustained,idx) = [UnitData(iun).Phase_spk(idx)];
% % % % % end
% % % % % 
% % % % % MeanPhase_Gap = nan(length(idxGap),5);
% % % % % for iun = idxGap'
% % % % %     idx = UnitData(iun).iSync;              
% % % % %     MeanPhase_Gap(iun==idxGap,idx) = [UnitData(iun).Phase_spk(idx)];
% % % % % end
% % % % % 
% % % % % 
% % % % % hfmp = figure;
% % % % % set(hfmp,'Position',tallrect)
% % % % % 
% % % % % for ir = 1:5
% % % % %     
% % % % %     isp(ir)=subplot(5,1,ir); hold on
% % % % %     ih(1,ir)=histogram( mod(MeanPhase_Gap(:,ir),360), 0:36:360,'FaceColor',colGap,'FaceAlpha',alphaGap);
% % % % %     ih(2,ir)=histogram( mod(MeanPhase_Sustained(:,ir),360), 0:36:360,'FaceColor',colSustained,'FaceAlpha',alphaSustained);
% % % % %     ih(3,ir)=histogram( mod(MeanPhase_Sparse(:,ir),360), 0:36:360,'FaceColor',colSparse,'FaceAlpha',alphaSparse);
% % % % %     xlim([0 360])
% % % % %     hold off
% % % % %     title([num2str(AMrates(ir)) ' Hz'])
% % % % % end
% % % % % 
% % % % % xlabel('Mean Phase (deg)')
% % % % % ylabel('Number of units')
% % % % % suptitle(sprintf('ex units, Sync stimuli only\nmean phase response'))
% % % % % 
% % % % % ymaxval = 10;
% % % % % linkaxes(isp,'xy')
% % % % % ylim([0 ymaxval])
% % % % % 
% % % % % legend(ih(:,5),{'Gap' 'Sustained' 'Sparse'},'Location','best')
% % % % % 
% % % % % % print_eps_kp(gcf,fullfile(savedir,'MeanPhaseDist_Sync'))
% % % % % % print_svg_kp(gcf,fullfile(savedir,'MeanPhaseDist_Sync'))
% % % % % 
% % % % % 






%% Waveform shapes

hf_tmp=figure;
SpikeWidthData   =[];

for iUn = 1:numel(UnitData)
    
%     % Only use labeled units
%     if strcmp(UnitInfo.RespType{iUn},'na') || strcmp(UnitInfo.RespType{iUn},'merged')
%         continue
%     end
    
    % Get this unit's info
    subject = UnitData(iUn).Subject;
    sess_tmp = strsplit(UnitData(iUn).Session,'_');
    session = sess_tmp{1};
    channel = UnitData(iUn).Channel(1);
    clu     = UnitData(iUn).Clu(1);
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
%         filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        clear Spikes Clusters
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Extract waveform info from each type of Spike data format
    %  should be interpolated to 
    
    if exist('Spikes','var') && ~exist('Clusters','var')         % >>> UMS <<<
        
        % Create intepolated version of waveform
        spikes = Spikes.sorted(channel);
        waveform = median(spikes.waveforms(spikes.assigns==clu,:),1);
                
    elseif exist('Clusters','var') && ~exist('Spikes','var')      % >>> KS <<<
        
        waveform = Clusters([Clusters.clusterID]==clu).maxChTemp;
        
    end
    
    waveform = waveform./abs(min(waveform));
    m = 40;
    x = (1:length(waveform)) /Info.fs*1000;
    q = linspace(min(x),max(x),length(x)*m);
    waveform_interp = interp1(x,waveform,q,'spline');
    
    
    % Now, if this was a merged unit, load the other waveform and compare
    % to the first
    if numel(sess_tmp)>1
        keyboard
        
        session2 = sess_tmp{2};
        clu2     = UnitData(iUn).Clu(2);
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session2); load(fullfile(fn.processed,subject,filename));
        
        % Create intepolated version of waveform
        spikes = Spikes.sorted(channel);
        waveform2 = median(spikes.waveforms(spikes.assigns==clu2,:),1);
        waveform2 = waveform2./abs(min(waveform2));
        m = 40;
        x = (1:length(waveform2)) /Info.fs*1000;
        q = linspace(min(x),max(x),length(x)*m);
        waveform_interp2 = interp1(x,waveform2,q,'spline');
        
        %         hf_tmp2 = figure; hold on
        %         plot(q,waveform_interp,'r','LineWidth',2)
        %         plot(q,waveform_interp2,'k','LineWidth',2)
        %
        %         keyboard
        %         close(hf_tmp2);
        
        waveform        = mean([waveform; waveform2],1);
        waveform_interp = mean([waveform_interp; waveform_interp2],1);
    end
    
    % Find peak and trough
    [trough,i_trough] = min(waveform_interp);    
    [~,i_peak] = max(waveform_interp);
        
    % Time from trough to peak (ms)
    ms_TroughPeak = q(i_peak)-q(i_trough);
    
    % Width (ms) at half amplitude
    halfmax = trough/2;
    [~,i_halfmax] = findpeaks( (-1*double(abs(waveform_interp-halfmax))) ,'SortStr','descend');
    iHalfMax = sort(i_halfmax(1:2));
    ms_WidthHalfMax = diff(q(iHalfMax));
    
    % Time to peak
%     [~,~,~,p] = findpeaks(double(-waveform_interp),'MinPeakProminence',abs(trough)/2);
%     [~,iprom] = min(abs(waveform_interp(1:i_trough)-(trough+p)));
    
    [~,ip,~,p] = findpeaks(double(waveform_interp));
    [~,iprom] = min(i_trough-ip(ip<i_trough));
    iprom = ip(iprom);
    if isempty(iprom), iprom=1; end
    ms_TimetoTrough = q(i_trough)-q(iprom);
    
    
    % Save waveform and metrics
    UnitData(iUn).waveform = waveform;
    SpikeWidthData = [SpikeWidthData; [  ms_WidthHalfMax  ms_TimetoTrough  ms_TroughPeak ]];
    
    
    % Plot if desired
%     figure(hf_tmp);
%     plot(x,waveform,'k')
%     hold on
%     plot(q,waveform_interp,'r','LineWidth',2)
%     plot([q(i_trough) q(i_peak)],[trough trough],'g','LineWidth',2)
%     plot(q,(-1*double(abs(waveform_interp-halfmax))) ,'k')
%     plot(q(iHalfMax),[halfmax halfmax],'b','LineWidth',2)
%     plot([q(iprom) q(i_trough)],[waveform_interp(iprom) waveform_interp(iprom)],'c','LineWidth',2)
%     clf
    
    
    
    
end %iUn


% Store waveform data in UnitInfo
UnitInfo.WidthHalfMax  = SpikeWidthData(:,1);
UnitInfo.TimetoTrough  = SpikeWidthData(:,2);
UnitInfo.TroughPeak    = SpikeWidthData(:,3);


% Plot waveform data
%   COLOR BY BASELINE FR
figure;
hold on
plot(UnitInfo.WidthHalfMax, UnitInfo.TroughPeak,'o','Color','k','LineWidth',2,'MarkerSize',10)
% plot(UnitInfo.WidthHalfMax(idxGap), UnitInfo.TroughPeak(idxGap),'o','Color',colGap,'LineWidth',2,'MarkerSize',10)
% plot(UnitInfo.WidthHalfMax(idxSustained), UnitInfo.TroughPeak(idxSustained),'o','Color',colSustained,'LineWidth',2,'MarkerSize',10)
% plot(UnitInfo.WidthHalfMax(idxSparse), UnitInfo.TroughPeak(idxSparse),'o','Color',colSparse,'LineWidth',2,'MarkerSize',10)
axis square
xlabel('Width at half max (ms)')
ylabel('Time from trough to peak (ms)')

% print_eps_kp(gcf,fullfile(savedir,'WaveformShapeMetrics'))
% print_svg_kp(gcf,fullfile(savedir,'WaveformShapeMetrics'))



figure;
subplot(2,2,1); hold on
plot(repmat(x,numel(idxGap),1)', vertcat(UnitData(idxGap).waveform)', 'Color',colGap,'LineWidth',2)
plot(repmat(x,numel(idxSustained),1)', vertcat(UnitData(idxSustained).waveform)', 'Color',colSustained,'LineWidth',2)
plot(repmat(x,numel(idxSparse),1)', vertcat(UnitData(idxSparse).waveform)', 'Color',colSparse,'LineWidth',2)
xlim([min(x) max(x)])
% ylim([-6 2].*10^-4)

subplot(2,2,2)
plot(repmat(x,numel(idxGap),1)', vertcat(UnitData(idxGap).waveform)', 'Color',colGap,'LineWidth',2)
xlim([min(x) max(x)])
% ylim([-6 2].*10^-4)

subplot(2,2,3)
plot(repmat(x,numel(idxSparse),1)', vertcat(UnitData(idxSparse).waveform)', 'Color',colSparse,'LineWidth',2)
xlim([min(x) max(x)])
% ylim([-6 2].*10^-4)

subplot(2,2,4)
plot(repmat(x,numel(idxSustained),1)', vertcat(UnitData(idxSustained).waveform)', 'Color',colSustained,'LineWidth',2)
xlim([min(x) max(x)])
% ylim([-6 2].*10^-4)

% print_eps_kp(gcf,fullfile(savedir,'Waveforms'))
% print_svg_kp(gcf,fullfile(savedir,'Waveforms'))



%%  Re-save UnitData and UnitInfo, with new data

% save(fullfile(fn.processed,'Units'),'UnitInfo','UnitData','-v7.3');




end



