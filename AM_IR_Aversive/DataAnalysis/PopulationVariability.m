function PopulationVariability
% 

% rate coding
% Fig 1: variability of avg FR per trial
%  5 sp (stimuli) FR range per cell
%  1 sp FF scatter plot

% temporal coding
% Fig 2: variability of timing of firing per trial
%  5 sp (stimuli) horizontal scatter of mean phase each trial, row per sig. cell
%  1 sp overall mean phase scatter plot

% Fig 3: more temporal coding plots
%  5 sp (stimuli) binarized MPH per sig. cell
%  5 sp (stimuli) zscored MPH per sig. cell 


% []for sig. sync. cell/rates, binarize MPH around baseline rate 
%    (demonstrates timing of responses across population)
% []for sig. sync. cell/rates, plot zscored MPH
% []plot mean phase per trial 
%    (analogous illustration to FR range, but shows variance of temporal coding across trials)

global AMrates 

alfaVS = 0.001;
trMax  = 40;
        
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

zFR_vec  = zeros(numel(UnitData),500,5);

Mus_2=[];  Uns_2=[];
Mus_3=[];  Uns_3=[];
Mus_4=[];  Uns_4=[];
Mus_5=[];  Uns_5=[];
Mus_6=[];  Uns_6=[];


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
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - 
    % - - - -   Some calculations require raw spike data   - - - - - 
    % - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - 
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
    FRmeds(iUn,:) = median(FRtrials(:,2:end),1);
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - -
    % - - - - - -     Variability of temporal coding     - - - - - - 
    % - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - 
    % Fig 2: variability of timing of firing per trial
    %  5 sp (stimuli) horizontal scatter of mean phase each trial, row per sig. cell
    %  1 sp overall mean phase scatter plot
    
    for stid = find(UnitData(iUn).VSdata_spk(2,:)>13.1)
        
        thisRaster = vertcat(MPH(MPH.ThisStimID==stid,:).raster{:});
        
        % Only analyze trials with spikes
        trials = find(sum(thisRaster,2)>0)';
        
        if isempty(trials)
            continue
        end
        
        % Randomize bc limiting to 40 plotted
        trials = trials(randperm(length(trials)));
        
        % Collect mean phase for each trial
        for it = trials(1:min(length(trials),40))
            
            Spktime=[];
            Spktime = find(thisRaster(it,:));
            
            mu = meanphase(Spktime,size(thisRaster,2));
            
            eval(sprintf('Mus_%i = [Mus_%i mu];',stid,stid))
            eval(sprintf('Uns_%i = [Uns_%i iUn];',stid,stid))
            
        end
        
        zFR_vec(iUn,1:size(thisRaster,2),stid-1) = mean(vertcat(MPH(MPH.ThisStimID==stid,:).zFR{:}),1);
        
        
    end %only significant sync
    
    
end %iUn



%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
twothirds   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)];
widescreen  = [1 scrsz(4) scrsz(3) scrsz(4)/3*2];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

ycountmax = 10*ceil(numel(UnitData)/2/10);

%%

%~~~~~~~~~~~~
%   Fig 1
%~~~~~~~~~~~~

hf1 = figure;
set(gcf,'Position',fullscreen)
hold on

% FR range
% [~,iu]=sort(FRmaxs(:,1)-FRmins(:,1));
[~,iu]=sort([UnitData.BaseFR]);
for ir=1:5
    
    subplot(3,5,ir+[0 5]);
    plot([FRmins(iu,ir) FRmaxs(iu,ir)]',[1:iUn; 1:iUn],'-k') %1+0.8.*[1:iUn; 1:iUn]./iUn-0.4
    hold on
    plot(FRmeds(iu,ir),1:iUn,'o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor','b')
    
    set(gca,'ytick',[],'xtick',0:20:80,'xticklabel',0:20:80,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    ylim([0 size(FRmeds,1)+1])
    xlim([0 80])
    if ir==1
        ylabel('Unit')
        xlabel('FR min to max')
    end
    title([num2str(AMrates(ir)) ' Hz'])
    box on
    axis fill
end


% FF: avg per stimulus
FFplot = nan(1,5);
for ist = 1:length(FFplot)
    FFplot(ist) = median(FFmeans(~isnan(FFmeans(:,ist+1)),ist+1));
end
maxFF=max(max(FFmeans(:,2:end)));

subplot(3,5,14:15);

plotSpread(FFmeans(:,2:end),'distributionColors','k')
plot(1:5,FFplot,'.-','MarkerSize',40,'LineWidth',2,'Color',[0.9 0.8 0])
set(gca,'xtick',1:5,'xticklabel',AMrates,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
xlim([-1 6])
ylim([0 10])
xlabel('AM rate (Hz)')
ylabel('Mean FF')
axis fill
box on



%~~~~~~~~~~~~
%   Fig 2
%~~~~~~~~~~~~

hf2 = figure;
set(gcf,'Position',fullscreen)
hold on

% Mean Phases

for ir=1:5
    
    subplot(2,5,ir+[0 5]);
    
    eval(sprintf('xdata = Mus_%i;',ir+1))
    eval(sprintf('ydata = Uns_%i;',ir+1))
    
    plot(xdata,ydata,'o','Color',[0 0 0.3],'MarkerSize',2)
    
    xlim([0 2*pi])
    ylim([0 iUn+1])
    set(gca,'ytick',[],'xtick',linspace(0,2*pi,5),'xticklabel',{'0' '' 'pi' '' '2pi'},'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    if ir==1
        ylabel('Unit')
        xlabel('Mean phase, each trial')
    end
    title([num2str(AMrates(ir)) ' Hz'])
    box on
    axis fill
    
end
    


%~~~~~~~~~~~~
%   Fig 3
%~~~~~~~~~~~~

hf3 = figure;
set(gcf,'Position',fullscreen)
hold on

% Mean Phases

for ir=1:5
    
    subplot(2,5,ir+[0 5]);
    
    zdata = zFR_vec(:,1:ceil(1000/AMrates(ir)),ir);
    
    imagesc(zdata)
    caxis([-1 3])
    cmocean('balance','pivot',0)
        
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0.5 iUn+0.5])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    if ir==1
        ylabel('Unit')
        xlabel('Time in MPH')
    end
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    
end





keyboard


% Save figures

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf1,fullfile(savedir,'VariabilitySummary_rate'))
print_eps_kp(hf2,fullfile(savedir,'VariabilitySummary_meanphase'))
print_eps_kp(hf3,fullfile(savedir,'VariabilitySummary_zscMPH'))


end





