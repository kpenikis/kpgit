function [UnitInfo,UnitData] = assessWaveformShapes(Units_fn, FigSavename, makePlot, UnitInfo, UnitData)
%
%  getWaveformShapes(UnitData,UnitInfo)
%    Calculates the spike width at half max and other parameters related to
%    unit type. Plot the population. Saves values into UnitInfo struct.
%  
%    Unit/response type not yet defined here. 
%
%  KP, 2019-02
%


% soon update to skip old iUn 



%% Load Unit data files

fn = set_paths_directories;

if nargin<5 && ~exist('UnitData','var')
        
    q = load(fullfile(fn.processed,[Units_fn '.mat']));
    UnitData = q.UnitData;
    UnitInfo = q.UnitInfo;
    clear q
    
    % [~,UnitData] = identifyResponsiveUnits(UnitData);
end
if isempty(FigSavename)
    FigSavename = [Units_fn '_WaveformShapes'];
end


%% Prepare figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
colormap('winter')

scrsz = get(0,'ScreenSize');
tallrect = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
widerect = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%% FR

AllStimResponses = cellfun(@(xxx) mean(xxx,1,'omitnan'), {UnitData.FR_raw_tr},'UniformOutput',false);
AllStimResponses = vertcat(AllStimResponses{:});

AllBaselineFRs = [UnitData.BaseFR]';

if ~isfield(UnitInfo,'TroughPeak')
    theseUnits = 1:size(UnitInfo,1);
else
    theseUnits = find(isnan([UnitInfo.TroughPeak]))';
end


%% Go through Units and get waveform shapes

SpikeWidthData=[];

for iUn = theseUnits
    
    % Get this unit's info
    subject = UnitData(iUn).Subject;
%     sess_tmp = strsplit(UnitData(iUn).Session,'_');
%     session = sess_tmp{1};
    session = UnitData(iUn).Session;
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
    
    % Find peak and trough
    [trough,i_trough] = min(waveform_interp);    
    [~,i_peak] = max(waveform_interp(i_trough:end));
    i_peak = i_peak + i_trough - 1;
        
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
    
    if ms_TroughPeak<0
        keyboard
    end
    
%     % Plot if desired
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


%%  Re-save UnitData and UnitInfo, with new data

UnitInfo.WidthHalfMax  = SpikeWidthData(:,1);
UnitInfo.TimetoTrough  = SpikeWidthData(:,2);
UnitInfo.TroughPeak    = SpikeWidthData(:,3);

% save(fullfile(fn.processed,Units_fn),'UnitInfo','UnitData','-v7.3');



%% Plot waveform distribution

if makePlot
    
    % Width x PeakTrough
% %     hf1 = figure;
% %     hold on
% %     scatter(UnitInfo.WidthHalfMax', UnitInfo.TroughPeak',200,log(ceil(AllBaselineFRs)'),'filled','Linewidth',3)
% %     % scatter(UnitInfo.WidthHalfMax', UnitInfo.TroughPeak', 300,'k')
% %     axis square
% %     set(gca,'xlim',[0.1 0.35],'ylim',[0.1 0.8])
% %     
% %     xlabel('Width at half max (ms)')
% %     ylabel('Time from trough to peak (ms)')
% %     title('Spike shape distribution, colored by baseline FR')
% %     text(0.27,0.15,sprintf('N = %i units\nfrom %i sessions',size(UnitInfo,1),numel(unique(UnitInfo.Session))))
% %     
% %     % colormap('winter')
% %     colormap(cmocean('-turbid'))
% %     cb=colorbar; pause(1)
% %     cb.Ruler.Scale = 'linear';
% %     cb.Ruler.TickValues = log(fliplr(ceil(max(AllBaselineFRs)):-5:min(AllBaselineFRs)));
% %     cb.Ruler.TickLabels = fliplr(ceil(max(AllBaselineFRs)):-5:min(AllBaselineFRs));
% %     cb.Label.String = 'Baseline FR (sp/s)';
% %     
% %     savedir = fullfile(fn.figs,'WaveformShapes');
% %     if ~exist(savedir,'dir')
% %         mkdir(savedir)
% %     end
% %     print_eps_kp(hf1,fullfile(savedir,FigSavename))
    
    
    % Baseline FR x PeakTrough
    
    hf2 = figure;
    hold on
    plot(AllBaselineFRs(AllBaselineFRs>0.1)', UnitInfo.TroughPeak(AllBaselineFRs>0.1)','.','Color',0.12*[1 1 1],'MarkerSize',20)
    plot(0.1*ones(sum(AllBaselineFRs<=0.1),1),UnitInfo.TroughPeak(AllBaselineFRs<=0.1)','ok')
    axis square
    set(gca,'xlim',[0.1 100],'ylim',[0.1 0.9],'xscale','log')
    
    xlabel('Spontaneous FR (sp/s)')
    ylabel('Time from trough to peak (ms)')
    title('Spike shape by baseline FR')
    text(0.27,0.15,sprintf('N = %i units\nfrom %i sessions',size(UnitInfo,1),numel(unique(UnitInfo.Session))))
    
    savedir = fullfile(fn.figs,'WaveformShapes');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(hf2,fullfile(savedir,[FigSavename '_bySpont']))
    
end


end



