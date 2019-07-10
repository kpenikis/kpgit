function PopulationTuning
% 

keyboard
% add subplot of zScore FRs, cdf/histograms for each stimulus

alfa = 0.05;
AMrates = [2 4 8 16 32];

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

%%

VS_sig = zeros(numel(UnitData),numel(AMrates));
BMF_FR = zeros(numel(UnitData),1);

for iUn = 1:numel(UnitData)
    
    % Percent of units synchronized by AMrate
    VS_sig(iUn,:) = UnitData(iUn).VSdata_spk(3,2:6)<alfa;
    
    % Distribution of BMF-FR and % of total units
    if ~isempty(UnitData(iUn).iBMF_FR)
        BMF_FR(iUn) = UnitData(iUn).iBMF_FR;
    end
    
end %iUn


%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',20)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
widerect   = [1 scrsz(4) scrsz(3)/3 scrsz(4)/4];


hf = figure;
set(gcf,'Position',widerect)

subplot(1,2,1)
plot(1:5,sum(VS_sig,1)/size(VS_sig,1),'.-k','MarkerSize',30,'LineWidth',2)
set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','Color','none')
xlim([0.5 5.5])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Proportion sig sync')

subplot(1,2,2)
histogram(BMF_FR,'Normalization','probability','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','Color','none','ytick',0:0.1:0.5)
xlim([-1 6])
ylim([0 0.5])
xlabel('BMF (calc via FR)')
ylabel('Proportion of units')


% Save figures

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,'TuningSummary'))


end







