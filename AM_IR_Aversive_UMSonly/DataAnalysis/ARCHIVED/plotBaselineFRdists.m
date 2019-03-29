function plotBaselineFRdists()
%
%  plotBaselineFRdists()
%    Goes through each unit in Resp struct and gets firing rate during
%    silent period at beginning of recording. Plots distribution of these
%    baseline FRs across the population, separating units labeled Sparse
%    from those labeled High FR.
%
%  KP, 2018-05
%


% Prepare figures
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
   

% Preallocate
BaselineFRs = [];
SpikeWidths = [];
AllWaves_Sparse = [];
AllWaves_HighFR = [];


% Load Resp data table 
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'RespStruct_allSU'));
Resp = q.Resp;

% Also load Redundant unit look up table
RedundantUnitLUT = readtable(fullfile(fn.processed,'RedundantUnitLUT'));


yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {Resp.FR_raw_tr},'UniformOutput',false);
Resp_meanFRs = mean(vertcat(yyy{:}),2,'omitnan');
[~,isrt] = sort(Resp_meanFRs);
Resp_meanFRs = Resp_meanFRs(isrt);
% Resp_sort = Resp(isrt);
Resp_sort = Resp;

theseUnits = 1:numel(Resp_sort);

% % Identify units labeled AM responsive
% [sigUnits,Resp_sort] = identifyResponsiveUnits(Resp_sort);


for iUn = theseUnits
    
    % Get this unit's info
    subject = Resp_sort(iUn).Subject;
    session = Resp_sort(iUn).Session;
    channel = Resp_sort(iUn).Channel;
    clu     = Resp_sort(iUn).Clu;
    
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    if (iUn>1 && ~( strcmp(subject,Resp_sort(iUn-1).Subject) && strcmp(session,Resp_sort(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,Resp_sort(iUn-1).Subject) && strcmp(session,Resp_sort(iUn-1).Session) && channel==Resp_sort(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes
    spikes = Spikes.sorted(channel);
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    
        
    % Skip if recorded more than once
%     if any( sum( [strcmp(RedundantUnitLUT.Subject,subject), strcmp(RedundantUnitLUT.Session,session),...
%              [RedundantUnitLUT.Channel]==channel, [RedundantUnitLUT.Clu]==clu ] , 2 ) == 4 )
%         continue
%     end
    
    
    %% Spike Width
    
    % Create intepolated version of waveform
    waveform = median(spikes.waveforms(spikes.assigns==clu,:),1);
    m = 40;
    x = (1:length(waveform)) /Info.fs*1000;
    q = linspace(min(x),max(x),length(x)*m);
    waveform_interp = interp1(x,waveform,q,'spline');
    
    % Save waveform for overlay plot
    if Resp_sort(iUn).Sparse==1
        AllWaves_Sparse = [AllWaves_Sparse; waveform_interp];
    else
        AllWaves_HighFR = [AllWaves_HighFR; waveform_interp];
    end
    
    % Find peak and trough
    [trough,i_trough] = min(waveform_interp);    
    [~,i_peak] = max(waveform_interp);
        
    % Time from trough to peak (ms)
    ms_TroughPeak = q(i_peak)-q(i_trough);
    
    % Width (ms) at half amplitude
    halfmax = trough/2;
    [peaks,i_halfmax] = findpeaks( (-1*double(abs(waveform_interp-halfmax))) );
    
    [peaks,ipeaks] = sort(peaks);
    peaks = peaks(end-1:end);
    ipeaks = ipeaks(end-1:end);
    
    iHalfMax = i_halfmax(sort(ipeaks));
    ms_WidthHalfMax = diff(q(iHalfMax));
    
    % Add to population data
    SpikeWidths = [SpikeWidths; [ ms_TroughPeak  ms_WidthHalfMax  Resp_sort(iUn).Sparse] ];
    
    
    % Plot if desired
%     figure(hf_tmp);
%     plot(x,waveform,'k')
%     hold on
%     plot(q,waveform_interp,'r','LineWidth',2)
%     plot([q(i_trough) q(i_peak)],[trough trough],'g','LineWidth',2)
%     plot(q,(-1*double(abs(waveform_interp-halfmax))) ,'k')
%     plot(q(iHalfMax),[halfmax halfmax],'b','LineWidth',2)
%     clf
%     
    
    %% Baseline FR
    
    % Convert FR to z-score
    bs_smth = 20;
    [~,~,Stream_Spikes,~] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    % Check that silence is first entry in TrialData
    if TrialData.trID(1)~=0
        keyboard
    end
    
    % Get onset and offset of silence
    t1 = TrialData.onset(1);
    t2 = TrialData.offset(1);
    
    % Calculate average spikes per second
    this_baseFR = sum(Stream_Spikes(t1+1:t2))/(t2-t1)*1000;
    
    % Add to population data
    BaselineFRs = [BaselineFRs; [ this_baseFR Resp_sort(iUn).Sparse] ];
    
    
end %iUn

keyboard 

%% Create figures

iHFR = SpikeWidths(:,3)==0;
iSps = SpikeWidths(:,3)==1;

% Baseline FR
hf1=figure;
hold on
ih(1)=histogram(BaselineFRs(iHFR,1), 0:5:40, 'FaceColor','k');
ih(2)=histogram(BaselineFRs(iSps,1), 0:5:40, 'FaceColor','b');
xlabel('FR (spk/s) during silence')
ylabel('Number of units')
legend(ih,{'High FR' 'Sparse'})
hold off

savedir = fullfile(fn.processed,'Rasters_SortFR');
% print_eps_kp(hf1,fullfile(savedir,'BaselineFR_Distributions_slim'))
% print_svg_kp(hf1,fullfile(savedir,'BaselineFR_Distributions_slim'))


% Spike Width  X  Baseline FR
hf=figure;
hold on
plot(SpikeWidths(iHFR,2),BaselineFRs(iHFR,1),'ok','LineWidth',2,'MarkerSize',10)
plot(SpikeWidths(iSps,2),BaselineFRs(iSps,1),'ob','LineWidth',2,'MarkerSize',10)
set(gca,'xscale','log','yscale','log','ytick',[0.1 1 2 5 10 20 40])
xlabel('Width at Half Max (ms)')
ylabel('Baseline FR (spk/s) during silence')
axis square

% print_eps_kp(hf,fullfile(savedir,'BaselineFR_SpikeWidth'))
% print_svg_kp(hf,fullfile(savedir,'BaselineFR_SpikeWidth'))


% Spike Width  X  Time to peak
hf2=figure;
hold on
plot(SpikeWidths(iHFR,2),SpikeWidths(iHFR,1),'ok','LineWidth',2)
plot(SpikeWidths(iSps,2),SpikeWidths(iSps,1),'ob','LineWidth',2)
xlim([0.1 0.3])
axis square
xlabel('Width at Half Max (ms)')
ylabel('Time from Trough to Peak (ms)')
title('Waveform metrics for all SU')

savedir = fullfile(fn.processed,'Rasters_SortFR');
print_eps_kp(hf2,fullfile(savedir,'SpikeWidth_Distributions_slim'))
print_svg_kp(hf2,fullfile(savedir,'SpikeWidth_Distributions_slim'))


% Now make figure overlaying waveforms for each type
keyboard
SEM_HighFR = std(AllWaves_HighFR,1)/sqrt(size(AllWaves_HighFR,1));
SEM_Sparse = std(AllWaves_Sparse,1)/sqrt(size(AllWaves_Sparse,1));

hf3=figure;
hold on
% plot(repmat(q,size(AllWaves_HighFR,1),1)',AllWaves_HighFR','k','LineWidth',0.15)
% plot(repmat(q,size(AllWaves_Sparse,1),1)',AllWaves_Sparse','b','LineWidth',0.15)
patch([q fliplr(q)],[mean(AllWaves_HighFR,1)+SEM_HighFR fliplr(mean(AllWaves_HighFR,1)-SEM_HighFR)], 'k', 'EdgeColor','none', 'FaceAlpha',0.4)
patch([q fliplr(q)],[mean(AllWaves_Sparse,1)+SEM_Sparse fliplr(mean(AllWaves_Sparse,1)-SEM_Sparse)], 'b', 'EdgeColor','none', 'FaceAlpha',0.4)
xlim([min(q) max(q)])
xlabel('ms')
ylim([-10 4].*10^-4)
title('Mean of waveforms, by type')
% title('All waveforms overlayed')


savedir = fullfile(fn.processed,'Rasters_SortFR');

print_eps_kp(hf3,fullfile(savedir,'WaveformsAll_slim'))
print_svg_kp(hf3,fullfile(savedir,'WaveformsAll_slim'))

print_eps_kp(hf3,fullfile(savedir,'WaveformsMean_slim'))
print_svg_kp(hf3,fullfile(savedir,'WaveformsMean_slim'))




end