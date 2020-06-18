function [UnitInfo,UnitData] = plotWaveformShapes(Units_fn)
%
%  plotWaveformShapes(UnitData,UnitInfo)
%    Calculates the spike width at half max and other parameters related to
%    unit type. Plot the population. Saves values into UnitInfo struct.
%  
%    Unit/response type not yet defined here. 
%
%  KP, 2020-04 
%


%% Load Unit data files

fn = set_paths_directories;

q = load(fullfile(fn.processed,[Units_fn '.mat']));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

AllBaselineFRs = [UnitData.BaseFR]';


%% Prepare figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
colormap('winter')

scrsz = get(0,'ScreenSize');     %[left bottom width height]
widerect = [1 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3];

savedir = fullfile(fn.figs,'WaveformShapes');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

FigSavename = [Units_fn '_WaveformShapes'];


%% Plot waveform distributions
    
% Baseline FR x PeakTrough

hf2 = figure;
set(gcf,'Position',widerect)

subplot(1,2,1)
hold on
plot(AllBaselineFRs(AllBaselineFRs>0.1)', UnitInfo.TroughPeak(AllBaselineFRs>0.1)','.','Color',0.12*[1 1 1],'MarkerSize',20)
plot(0.1*ones(sum(AllBaselineFRs<=0.1),1),UnitInfo.TroughPeak(AllBaselineFRs<=0.1)','ok')
axis square
set(gca,'Color','none','xlim',[0.1 100],'ylim',[0.1 0.9],'xscale','log')

xlabel('Spontaneous FR (sp/s)')
ylabel('Time from trough to peak (ms)')
% title('Spike shape by baseline FR')
text(0.27,0.15,sprintf('N = %i units\nfrom %i sessions',size(UnitInfo,1),numel(unique(UnitInfo.Session))))


subplot(1,2,2)
histogram(UnitInfo.TroughPeak,0.1:0.03333:0.9)
% [caz,cel] = view;
set(gca,'Color','none','ylim',[0 size(UnitInfo,1)/3],'xlim',[0.1 0.9])
view(90,-90)

print_eps_kp(hf2,fullfile(savedir,[FigSavename '_bySpont']))



    
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
% %     print_eps_kp(hf1,fullfile(savedir,FigSavename))
    
    
end



