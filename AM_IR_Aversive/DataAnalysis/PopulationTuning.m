function PopulationTuning
% 


% add subplot of zScore FRs, cdf/histograms for each stimulus

alfa = 0.001;
AMrates = [2 4 8 16 32];

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
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load info struct
filename = sprintf( '%s_sess-%s_Info',UnitData(1).Subject,UnitData(1).Session); 
load(fullfile(fn.processed,UnitData(1).Subject,filename));


%%

VS_sig = zeros(numel(UnitData),numel(AMrates));
BMF_FR = zeros(numel(UnitData),1);
zFRs   = nan(numel(UnitData),7);
dFRs   = nan(numel(UnitData),7);

for iUn = 1:numel(UnitData)
    
    % Percent of units synchronized by AMrate
    VS_sig(iUn,:) = UnitData(iUn).VSdata_spk(3,2:6)<alfa;
    
    % Distribution of BMF-FR and % of total units
    if ~isempty(UnitData(iUn).iBMF_FR)
        BMF_FR(iUn) = UnitData(iUn).iBMF_FR;
    end
    
    % Avg FR per stimulus, z-scored
    zFRs(iUn,:) = UnitData(iUn).FR_nrm(2:end);
    
    dFRs(iUn,:) = mean(UnitData(iUn).FR_raw_tr(:,2:end),1,'omitnan') - UnitData(iUn).BaseFR;
    
end %iUn


%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
widerect   = [1 scrsz(4) scrsz(3) scrsz(4)/3];


hf = figure;
set(gcf,'Position',widerect)

subplot(1,5,1)
plot(1:5,sum(VS_sig,1)/size(VS_sig,1),'.-k','MarkerSize',50,'LineWidth',3)
set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','Color','none')
xlim([0.5 5.5])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Proportion sig sync')
axis square

subplot(1,5,2)
histogram(BMF_FR,'Normalization','probability','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','Color','none','ytick',0:0.1:0.5)
xlim([-1 6])
ylim([0 0.5])
xlabel('BMF (calc via FR)')
ylabel('Proportion of units')
axis square


subplot(1,5,3)
hold on
for ist = 1:size(zFRs,2)
    ih(ist) = histogram(zFRs(:,ist),-5:0.001:10,'Normalization','cdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
%     ih(ist) = plot(-5:0.001:10,foo,'Color',colors(ist+1,:),'LineWidth',2);
end
axis square
xlim([-1 1])
xlabel('FR (z-score)')
ylabel('Cumulative probability')
legend(ih,Info.stim_ID_key(2:end),'FontSize',9,'Location','southeast')
legend('boxoff')
grid on

subplot(1,5,4)
set(gca,'xscale','log')
hold on
ih0=histogram([UnitData.BaseFR],logspace(log10(0.1),log10(100),100),'Normalization','cdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',4);
for ist = 1:size(zFRs,2)
    ih(ist) = histogram(dFRs(:,ist)+[UnitData.BaseFR]',logspace(log10(0.1),log10(100),100),'Normalization','cdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
end
axis square
xlim([0.1 100])
ylim([0 1])
xlabel('FR (Hz)')
ylabel('Cumulative probability')
legend(ih0,{'Silence'},'FontSize',9,'Location','southeast')
legend('boxoff')
grid on


subplot(1,5,5)
hold on
for ist = 1:size(zFRs,2)
    ih(ist) = histogram(dFRs(:,ist),-50:0.1:100,'Normalization','cdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
end
axis square
xlim(20*[-1 1])
ylim([0 1])
xlabel('Sound driven FR (Resp-BaseFR)')
ylabel('Cumulative probability')
legend(ih,Info.stim_ID_key(2:end),'FontSize',9,'Location','southeast')
legend('boxoff')
grid on



% Save figures

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,'TuningSummary'))


end







