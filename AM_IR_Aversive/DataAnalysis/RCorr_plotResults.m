function RCorr_plotResults


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

ntTemp = 1;
% q = load(fullfile(fn.figs,'RCorr','exclOnset','PCMat_10trTemp'));
q = load(fullfile(fn.figs,'RCorr','exclOnset',sprintf('PCMat_%itrT',ntTemp)));
PCMat = q.PCMat;
clear q


% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
twothirds   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)];


Thresholds = 0.6:-0.1:0.5;

for iTh = 1:numel(Thresholds)
    AccuracyThreshold = Thresholds(iTh);
    
    nStim_Classified = nan(size(PCMat,3),1);
    
    for iUn = 1:size(PCMat,3)
        nStim_Classified(iUn) = sum(diag(PCMat(:,:,iUn))>AccuracyThreshold);
    end
    
    hf(iTh)=figure;
    histogram(nStim_Classified,'FaceColor',[1 1 1].*0.7)
    hold on
    text(-0.2,18,num2str(sum(nStim_Classified==0)))
    ylim([0 20])
    xlim([-1 9])
    set(gca,'xtick',0:8,'tickdir','out')
    title(['Number of stim classfied above ' num2str(AccuracyThreshold) ' PC'])
    
    print_eps_kp(hf(iTh), fullfile(fn.figs,'RCorr','exclOnset',sprintf('nStimHist_%itrT_%i',ntTemp,10*AccuracyThreshold)))
    
end

end
