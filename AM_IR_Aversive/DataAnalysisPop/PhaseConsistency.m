% for each cell
% mean phase across rates (indicate if significant or not)
% phase plot for each cell? 
% single measure to pull out: 
%  -proportion of significant rates within 90 degrees
%  -weighted vector

function PhaseConsistency
%
% PhaseConsistency
%
%  Code orginated from savePopMPHdata. To send to Yulia.
%

global AMrates

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];


%% Plot raw data: concentric sig response phases 

PhaseData = nan(numel(UnitData),3);

plotcolors = cmocean('phase',361);

figure;
% set(gcf,'Position',fullscreen)
hp=polaraxes;
set(hp,'Color','none')
hold on

for iUn = 1:numel(UnitData)
    
    % Find significantly synchronized phases
    if sum(UnitData(iUn).VSdata_spk(2,2:4)>13.1)<3
        continue
    end
    
    PhaseData(iUn,:) = UnitData(iUn).Phase_spk(2:4);
    
    % Plot
    polarplot([0 fliplr(deg2rad(PhaseData(iUn,:)))],[0 1 2 3],...
        '-','Color',plotcolors(round(PhaseData(iUn,1))+1,:),...
        'LineWidth',2)
       
end %iUn

title('Phases of responses to 2,4,8 Hz, for each unit with all sig sync')
rlim([0 3])
set(gca,'rtick',[])
hp.ThetaZeroLocation = 'bottom';
hp.ThetaDir = 'clockwise';

% Save figure
savedir = fullfile(fn.figs,'PhaseShifts');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

print(fullfile(savedir,'SigPhases_2-4-8'),'-dpdf','-bestfit')


%% Plot raw data: sig response phases 

% PhaseData = nan(numel(UnitData),5);
% 
% figure;
% set(gcf,'Position',fullscreen)
% hp=polaraxes;
% set(hp,'Color','none')
% hold on
% 
% for iUn = 1:numel(UnitData)
%     
%     PhaseData(iUn,:) = UnitData(iUn).Phase_spk(2:6);
%     
%     % Wrap latency around to next cycle where necessary   %%not needed for polar plot
%     %     resets = 1+find(diff(PhaseData(iUn,:))<30);
%     %     for ii = 1:numel(resets)
%     %         PhaseData(iUn,resets(ii):end) = PhaseData(iUn,resets(ii):end) + 2*pi;
%     %     end
%     
%     % Find significantly synchronized phases
%     iSig = UnitData(iUn).VSdata_spk(2,2:6)>13.1;
%     
%     % Plot
% %     polarplot(deg2rad(PhaseData(iUn,iSig)),UnitData(iUn).VSdata_spk(1,find(iSig)+1),...
% %         '-','Color',[0.071 0.573 0.686],'LineWidth',1)
%     polarplot(deg2rad(PhaseData(iUn,iSig)),ones(size(PhaseData(iUn,iSig))),...
%         '-','LineWidth',1)
%        
% end %iUn
% 
% title('Phases of sig sync responses, connected, for each unit')
% rlim([0 1])
% hp.ThetaZeroLocation = 'bottom';
% hp.ThetaDir = 'clockwise';
% 
% % Save figure
% savedir = fullfile(fn.figs,'PhaseShifts');
% if ~exist(savedir,'dir')
%     mkdir(savedir)
% end
% 
% print(fullfile(savedir,'AllSigPhaseShifts'),'-dpdf','-bestfit')



%% Histogram of PHASES for each rate

% 3rd dim: [mu VS RS]
PhaseData = nan(numel(UnitData),5,3);

for iUn = 1:numel(UnitData)
    PhaseData(iUn,:,1) = UnitData(iUn).Phase_spk(2:6);
    PhaseData(iUn,:,2) = UnitData(iUn).VSdata_spk(1,2:6);
    PhaseData(iUn,:,3) = UnitData(iUn).VSdata_spk(2,2:6);
end %iUn


% For each doubling of AM rate, with both sig sync, get distribution of
% phase shifts

Theta = linspace(0,2*pi,16);
plotTheta = diff(Theta)./2 + Theta(1:end-1);

iSig = PhaseData(:,:,3)>13.1;

hf=figure; 
set(hf,'Position',fullscreen)

for ir = 1:size(PhaseData,2)
    
    Phases = PhaseData( sum(iSig(:,ir),2)==1 ,ir,1);
    HistData = histcounts(deg2rad(Phases),Theta);
    
    % Plot result
    subplot(2,3,ir,polaraxes); 
    polarplot([plotTheta plotTheta(1)],[HistData HistData(1)],'Color',[0.173 0.171 0.686],'LineWidth',2)
    set(gca,'Color','none','RAxisLocation',130)
    set(gca,'ThetaZeroLocation','bottom',...
        'ThetaDir','clockwise');
    title(sprintf('Distribution of phases at %i Hz',AMrates(ir)))
    set(gca,'rlim',[0 35])
    
end


% Save figure
set(gcf,'PaperOrientation','landscape')
print(fullfile(savedir,'PhaseDistributions'),'-dpdf','-bestfit')



%% Histogram of phase SHIFTS with doubling rate

% 3rd dim: [mu VS RS]
PhaseData = nan(numel(UnitData),5,3);

for iUn = 1:numel(UnitData)
    PhaseData(iUn,:,1) = UnitData(iUn).Phase_spk(2:6);
    PhaseData(iUn,:,2) = UnitData(iUn).VSdata_spk(1,2:6);
    PhaseData(iUn,:,3) = UnitData(iUn).VSdata_spk(2,2:6);
end %iUn


% For each doubling of AM rate, with both sig sync, get distribution of
% phase shifts

Theta = linspace(-2*pi,2*pi,32);
plotTheta = diff(Theta)./2 + Theta(1:end-1);

iSig = PhaseData(:,:,3)>13.1;

hf=figure; 
set(hf,'Position',fullscreen)

for ir = 1:size(PhaseData,2)-1
    
    PhaseShifts = diff( PhaseData( sum(iSig(:,ir:ir+1),2)==2 ,ir:ir+1,1), 1,2);
    HistData = histcounts(deg2rad(PhaseShifts),Theta);
    
    % Plot result
    subplot(2,2,ir,polaraxes);
    polarplot(plotTheta,HistData,'Color',[0.071 0.686 0.573],'LineWidth',2)
    set(gca,'Color','none','RAxisLocation',130)
    set(gca,'ThetaZeroLocation','bottom',...
        'ThetaDir','clockwise');
    title(sprintf('%i Hz phase - %i Hz phase',AMrates(ir+1),AMrates(ir)))
    set(gca,'rlim',[0 60])
end


% Save figure
set(gcf,'PaperOrientation','landscape')
print(fullfile(savedir,'PhaseSHIFTDistributions'),'-dpdf','-bestfit')



%% Correlation of phase SHIFT with doubling rate

% 3rd dim: [mu VS RS]
PhaseData = nan(numel(UnitData),5,3);

for iUn = 1:numel(UnitData)
    PhaseData(iUn,:,1) = UnitData(iUn).Phase_spk(2:6);
    PhaseData(iUn,:,2) = UnitData(iUn).VSdata_spk(1,2:6);
    PhaseData(iUn,:,3) = UnitData(iUn).VSdata_spk(2,2:6);
end %iUn


% For each doubling of AM rate, with both sig sync, get distribution of
% phase shifts

Theta = linspace(-2*pi,2*pi,32);
plotTheta = diff(Theta)./2 + Theta(1:end-1);

iSig = PhaseData(:,:,3)>13.1;

hf=figure; 
set(hf,'Position',fullscreen)

for ir = 1:size(PhaseData,2)-1
    
    % Plot result
    subplot(2,2,ir);
    plot([0 360],[0 360],'Color',0.7*[1 1 1],'LineWidth',2)
    hold on
    plot(PhaseData( sum(iSig(:,ir:ir+1),2)==2 ,ir,1),PhaseData( sum(iSig(:,ir:ir+1),2)==2 ,ir+1,1),'.','Color',[0.071 0.686 0.573],'MarkerSize',20)
    axis square
    xlabel(sprintf('%i Hz phase',AMrates(ir)))
    ylabel(sprintf('%i Hz phase',AMrates(ir+1)))
    xlim([0 360])
    ylim([0 360])
    
end


% Save figure
set(gcf,'PaperOrientation','landscape')
print(fullfile(savedir,'PhaseCorrelations'),'-dpdf','-bestfit')




end




