function PopulationTuning
% 


% [ ] proportion units significantly above silence baseline 
% [ ] distribution of BMF-VS
% [ ] cdf to pdf, don?t plot z-score
% [ ] trial variability plots 
% [ ] bar: mean FF as function of stimulus (last 750 ms)
% [ ] scatter: FF as a function of FR  (last 750 ms)
% [ ] median CorrScore as a function of stimulus (last 750 ms)
% [ ] mean percentage ?miss? trials as a function of stimulus (last 750 ms)


% add subplot of zScore FRs, cdf/histograms for each stimulus

alfaVS = 0.001;
alfaFR = 0.05;
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
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%%

VS_sig   = zeros(numel(UnitData),numel(AMrates));
BMF_VS   = zeros(numel(UnitData),1);

FR_sig05 = zeros(numel(UnitData),numel(AMrates));
FR_sig01 = zeros(numel(UnitData),numel(AMrates));
FR_dp    = zeros(numel(UnitData),numel(AMrates));
BMF_FR   = zeros(numel(UnitData),1);

FRmins   = nan(numel(UnitData),numel(AMrates));
FRmaxs   = nan(numel(UnitData),numel(AMrates));

FFmeans  = nan(numel(UnitData),numel(AMrates)+1);

zFRs     = nan(numel(UnitData),7);
dFRs     = nan(numel(UnitData),7);

for iUn = 1:numel(UnitData)
    
    % Percent of units synchronized by AMrate
    VS_sig(iUn,:) = UnitData(iUn).VSdata_spk(3,2:6)<alfaVS;
    
    % Distribution of BMF-VS and N ns
    if ~isempty(UnitData(iUn).iBMF_VS)
        BMF_VS(iUn) = UnitData(iUn).iBMF_VS;
    end
    
    % Percent of units FR > silence per stimulus, dprime
    FR_dp(iUn,:) = abs(UnitData(iUn).dp_mat(2:6,2))';
    
    % Distribution of BMF-FR and N ns
    if ~isempty(UnitData(iUn).iBMF_FR)
        BMF_FR(iUn) = UnitData(iUn).iBMF_FR;
    end    
    
    % Avg FR per stimulus, z-scored
    zFRs(iUn,:) = UnitData(iUn).FR_nrm(2:end);
    
    dFRs(iUn,:) = mean(UnitData(iUn).FR_raw_tr(:,2:end),1,'omitnan') - UnitData(iUn).BaseFR;
    
    
    % - - - -  The following calculations require raw spike data - - - - - 
    get_trial_data_posthoc
    
    
    % Percent of units FR significantly > silence
    
    pvals = nan(1,size(FRtrials,2)-1);
    for ist = 1:length(pvals)
        pvals(ist) = ranksum(FRtrials(:,1),FRtrials(:,ist+1));
    end
    
    [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,0.05);
    FR_sig05(iUn,:) = iIR_sig_bonferroni;
    
    [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,0.01);
    FR_sig01(iUn,:) = iIR_sig_bonferroni;
    
    
    % Fano factor for each stimulus (last 750 ms)
    FFmeans(iUn,:) = var(FRtrials,1)./mean(FRtrials,1);
    
    % Range of FR across trials
    FRmins(iUn,:) = min(FRtrials(:,2:end),[],1); %min(UnitData(iUn).FR_raw_tr(:,2:6))
    FRmaxs(iUn,:) = max(FRtrials(:,2:end),[],1); %max(UnitData(iUn).FR_raw_tr(:,2:6))
    
end %iUn


%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
widerect    = [1 scrsz(4) scrsz(3) scrsz(4)/3];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

ycountmax = 10*ceil(numel(UnitData)/2/10);

hf = figure;
set(gcf,'Position',fullscreen)

% VS: Prop Sync
subplot(3,4,2)
plot(1:5,sum(VS_sig,1)/size(VS_sig,1),'.-r','MarkerSize',50,'LineWidth',4)
set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
xlim([-1 6])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Proportion sig sync')
axis square
box on

% VS: BMF distribution
subplot(3,4,6)
histogram(BMF_VS,'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax/2 ycountmax])
xlim([-1 6])
ylim([0 ycountmax])
xlabel('BMF (VS)')
ylabel('N units')
axis square
box on

% FR: Prop diff silence
subplot(3,4,3)
hold on
% dprime formula
plot(1:5,sum(FR_dp>1,1)/size(FR_dp,1),'-','LineWidth',4,'Color',[148 251 145]./255)
plot(1:5,sum(FR_dp>1.5,1)/size(FR_dp,1),'-','LineWidth',4,'Color',[84 191 97]./255)
plot(1:5,sum(FR_dp>2,1)/size(FR_dp,1),'-','LineWidth',4,'Color',[29 133 57]./255)
% ranksum test
plot(1:5,sum(FR_sig05,1)/size(FR_sig05,1),'-','LineWidth',4,'Color',[0.3 0.5 0.95])
plot(1:5,sum(FR_sig01,1)/size(FR_sig01,1),'-','LineWidth',4,'Color',[0.2 0.33 0.8])
set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
xlim([-1 6])
ylim([0 1])
xlabel('AM rate (Hz)')
ylabel('Proportion FR>silence')
axis square
box on

% VS: BMF distribution
subplot(3,4,7)
histogram(BMF_FR,'Normalization','count','FaceColor','k','EdgeColor','none','FaceAlpha',1);
set(gca,'xtick',0:5,'xticklabel',{'ns' AMrates},'tickdir','out','ticklength',[0.02 0.02],'Color','none','ytick',[0 ycountmax/2 ycountmax])
xlim([-1 6])
ylim([0 ycountmax])
xlabel('BMF (FR)')
ylabel('N units')
axis square
box on

% 
% % FF: avg per stimulus
% FFplot = nan(1,5);
% for ist = 1:length(FFplot)
%     FFplot(ist) = median(FFmeans(~isnan(FFmeans(:,ist+1)),ist+1));
% end
% maxFF=max(max(FFmeans(:,2:end)));
% 
% subplot(3,4,10)
% plotSpread(FFmeans(:,2:end),'distributionColors','k')
% plot(1:5,FFplot,'.-','MarkerSize',50,'LineWidth',4,'Color',[0.9 0.8 0])
% set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
% xlim([-1 6])
% ylim([0 10])
% xlabel('AM rate (Hz)')
% ylabel('Mean FF')
% axis square
% box on
% 
% 
% % FR range
% % [~,iu]=sort(FRmaxs(:,1)-FRmins(:,1));
% [~,iu]=sort([UnitData.BaseFR]);
% subplot(3,4,[1 5 9]);
% for ir=1:5
%     plot([FRmins(iu,ir) FRmaxs(iu,ir)]',ir+0.8.*[1:iUn; 1:iUn]./iUn-0.4,'-k')
%     hold on
% end
% % plot(iu,FFplot,'.-','MarkerSize',50,'LineWidth',4,'Color',[0.9 0.8 0])
% set(gca,'ytick',1:5,'yticklabel',AMrates,'xtick',0:20:80,'xticklabel',0:20:80,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
% ylim([0.5 5.5])
% xlim([0 80])
% ylabel('AM rate (Hz)')
% xlabel('FR min to max')
% box on
% axis fill


% Distribution of FRs per stimulus
subplot(3,4,4)
set(gca,'xscale','log')
hold on
for ist = 1:size(zFRs,2)
    ih(ist) = histogram(dFRs(:,ist)+[UnitData.BaseFR]',logspace(log10(0.1),log10(100),410),'Normalization','cdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
end
ih0=histogram([UnitData.BaseFR],logspace(log10(0.1),log10(100),410),'Normalization','cdf','DisplayStyle','stairs','EdgeColor','k','LineWidth',4);
axis square
set(gca,'Color','none')
xlim([0.1 100])
ylim([0 1])
xlabel('FR Resp (sp/s)')
ylabel('Cumulative probability')
legend(ih0,{'Silence'},'FontSize',9,'Location','southeast')
legend('boxoff')
grid on


% Distribution of baseline-subtracted FRs per stimulus
subplot(3,4,8)
% hold off; cla

window=gausswin(10);
window=window-min(window);
window=window/sum(window);

for ist = 1:size(dFRs,2)
%     ih(ist) = histogram(dFRs(:,ist),-50:0.1:100,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',colors(ist+1,:),'LineWidth',2);
    [values,edges] = histcounts(dFRs(:,ist),-31:0.5:31);
    sig=conv(values,window,'same');
    ih(ist) = plot(edges(2:end)-mode(diff(edges))/2,sig./sum(sig),'Color',colors(ist+1,:),'LineWidth',2);
    hold on
end
axis square
set(gca,'Color','none')
xlim(30*[-1 1])
ylim([0 0.1])
xlabel('Sound driven FR')
ylabel('Probability')
legend(ih,Info.stim_ID_key(2:end),'FontSize',9,'Location','northeast')
legend('boxoff')
grid on



% Save figures

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,'TuningSummary'))


end







