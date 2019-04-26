function playNoiseCorr(SUBJECT,SESSION)
% 
% Must be single session.
% Gathers raster data for all units and all stimuli. 
% Then look at them in a bunch of different ways.
% 
% KP, 2019-03
% 


close all
fn = set_paths_directories(SUBJECT,SESSION,1);

%%
% Get raster data
[UnitInfo, UnitData, Info, TrialData, Clusters, StimResp ] = collectRasterDataSession(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

Clu1 = 1230;
Clu2 = 68;
Un1  = [UnitInfo.Clu]==Clu1;
Un2  = [UnitInfo.Clu]==Clu2;
spiketimes1 = Clusters([Clusters.clusterID]==Clu1).spikeTimes;
spiketimes2 = Clusters([Clusters.clusterID]==Clu2).spikeTimes;

fsm = 10;
convWin = 500*fsm;
recLen = max([TrialData.offset(end) max(spiketimes1) max(spiketimes2)]);

sp1 = zeros(1,recLen*fsm);
sp1(round(spiketimes1*1000*fsm)) = 1;
sp1 = convolveGauss(sp1,convWin);

sp2 = zeros(1,recLen*fsm);
sp2(round(spiketimes2*1000*fsm)) = 1;
sp2 = convolveGauss(sp2,convWin);


winLen    = 1000*fsm;
winStp    =  1*fsm;
t0        = winLen/2;
winCents  = t0:winStp:(recLen*fsm-t0);
nWins     = numel(winCents);

runningRnoise = nan(1,nWins);
for iWin = 1:nWins
    t1 = 1 + winCents(iWin) - winLen/2;
    t2 = t1 + winLen-1;
    r = corrcoef(sp1(t1:t2),sp2(t1:t2));
    runningRnoise(1,iWin) = r(1,2);
end

runningRnoise = [nan(1,t0/fsm+1) runningRnoise nan(1,t0/fsm)];

% figure;
% plot([0 recLen*fsm],[0 0],'k')
% hold on
% plot( winCents, runningRnoise, 'b','LineWidth',2)
% xlim([0 recLen*fsm])


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
rng('shuffle')

scrsz = get(0,'ScreenSize');    %[left bottom width height]
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];
tallrect  = [1 scrsz(4) scrsz(3)/3 scrsz(4)];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937];


% Set colors
colors = [   0 200 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        
rasterdotsize  = 18;
psthlinewidth  = 4;


%% Now collect R noise values during each stimulus

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials] = get_clean_trials( TrialData,...
    unique([Info.artifact(UnitData(Un1).Channel).trials; Info.artifact(UnitData(Un2).Channel).trials]),...
    UnitData(Un1).spl,UnitData(Un1).lpn);
allStim = unique(TrialData.trID(all_TDidx))';

minTrs = min(Ntrials(~isnan(Ntrials)));

ymaxval = 0;

Rnoise_RUN_stim     = nan(1,numel(allStim));
Rnoise_RUN_stim_std = nan(1,numel(allStim));

Rnoise_stim       = nan(1,numel(allStim));

FR1_stim          = nan(1,numel(allStim));
FR2_stim          = nan(1,numel(allStim));

zFR1 = [];
zFR2 = [];

for ist = 1:numel(allStim)
    
    stid = allStim(ist);
    
    %% Collect trial indices and timestamps
    
    if stid==3 || stid==6
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
        % Find Pdc trials that follow same rate during ITI
        TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
        
        TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
        TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        
    else
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
        TDidx_iti = [];
    end
    
    
    % Get timestamps of onsets and offsets
    clear t2 t3 Duration t_win
    t2 = TrialData.onset(TDidx);
    t3 = TrialData.offset(TDidx);
    Duration = mode(diff([t2 t3],1,2));
    
    % Add ITI trials (shortened to match duration)
    if ~isempty(TDidx_iti)
        t2 = [t2; TrialData.onset(TDidx_iti)];
        TDidx = [TDidx; TDidx_iti];
    end
    
%     kt = randperm(length(t2));
%     kt = kt(1:minTrs);
    kt = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    
    %%
    % Set up figure
    hf(ist) = figure;
    set(gcf,'Position',tallrect.* [1 1 figwidthscales(stid) 1])
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
    hs(ist,2) = subplot(9,1,2:5);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[],...
        'ylim', [0 length(t2)+1])
    hs(ist,3) = subplot(9,1,6:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'xticklabel',[0 (figwidthscales(stid)*1000)])
    
    
    %%
    % Preallocate
    stim   = nan( numel(TDidx), Duration+1 );
    Rnoise = nan( numel(TDidx), Duration+1 );
    FR1 = [];
    FR2 = [];
    
    % Collect rms/noise for this stimulus/unit
    for it = 1:numel(TDidx)
        
        Rnoise(it,:) = ...
            runningRnoise( (t2(it) : t3(it)) );
        
        stim(it,:) = ...
            SoundStream(1, t2(it) : t3(it) )...
            ./ max(SoundStream(1, t2(it) : t3(it) ));
        
        FR1 = [FR1 sum(sp1((fsm*t2(it)+1):(fsm*t3(it))))];
        FR2 = [FR2 sum(sp2((fsm*t2(it)+1):(fsm*t3(it))))];
        
    end %it
    
    
    %% Add to plots
    
    figure(hf(ist)); hold on
    
    % Stimulus
    subplot(hs(ist,1)); hold on
    plot(0:Duration, mean(stim,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color',colors(stid,:))
    set(gca,'Color','none','xtick',[],'ytick',[])
    hold off
    
    % R noise -- raw trajectories
    subplot(hs(ist,2)); 
    plot(repmat(0:Duration,size(Rnoise,1),1)', Rnoise',...
        'LineWidth',1,'Color','k')
    hold on
    plot(0:Duration, mean(Rnoise,1,'omitnan'),...
        'LineWidth',2*psthlinewidth,'Color',colors(stid,:))
    set(gca,'Color','none','xtick',[],'tickdir','out')
    box off
    hold off
    
    %  R noise -- norm 
    subplot(hs(ist,3)); hold on
    plot(repmat(0:Duration,size(Rnoise,1),1)', Rnoise' - repmat(Rnoise(:,1)',size(Rnoise,2),1),...
        'LineWidth',1,'Color','k')
    plot(0:Duration, mean(Rnoise' - repmat(Rnoise(:,1)',size(Rnoise,2),1),2,'omitnan')',...
        'LineWidth',2*psthlinewidth,'Color',colors(stid,:))
    set(gca,'Color','none',...
        'tickdir','out')
    hold off
    
    ymaxval = max(ymaxval,max(mean(Rnoise,1,'omitnan')));
    
    suptitle(sprintf('%s\n%s %s\nNoise Correlations | Clus: %i & %i',...
        Info.stim_ID_key{stid}, SUBJECT,SESSION,Clu1,Clu2 ))
    
    % Save fig
    savedir = fullfile(fn.figs,'NoiseCorr','Running');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(hf(ist),fullfile(savedir,sprintf('%s_%s_%s_Clus_%i_%i_500msConv',SUBJECT,SESSION,Info.stim_ID_key{ist},Clu1,Clu2)))
    
    %%
    
    Rnoise_RUN_stim(ist)     = mean(mean(Rnoise,2,'omitnan'),'omitnan');
    Rnoise_RUN_stim_std(ist) = std(mean(Rnoise,2,'omitnan'),'omitnan');
    
    r = corrcoef(FR1,FR2);
    Rnoise_stim(ist)       = r(1,2);
    
    FR1_stim(ist) = mean(FR1);
    FR2_stim(ist) = mean(FR2);
    
    zFR1 = [zFR1 zscore(FR1)];
    zFR2 = [zFR2 zscore(FR2)];
    
end %ist


figure; hold on

% Signal correlation
r = corrcoef(FR1_stim,FR2_stim);
plot([0 ist+1],[r(1,2) r(1,2)],'--k','LineWidth',2)

% Noise correlation (all stim, method downer paper)
r = corrcoef(zFR1,zFR2);
plot([0 ist+1],[r(1,2) r(1,2)],'--m','LineWidth',2)

% Noise correlation (whole session)
minL = min([length(sp1) length(sp2)]);
r = corrcoef(sp1(TrialData.onset(1):minL),sp2(TrialData.onset(1):minL));
plot([0 ist+1],[r(1,2) r(1,2)],'--r','LineWidth',2)

% Noise correlation (silence)
r = corrcoef(sp1(TrialData.onset(1):TrialData.offset(1)),sp2(TrialData.onset(1):TrialData.offset(1)));
plot([0 ist+1],[r(1,2) r(1,2)],'--b','LineWidth',2)

% Each stimulus noise correlation - from running corr
for ist = 1:numel(allStim)
    plot([allStim(ist); allStim(ist)],Rnoise_RUN_stim(ist)+Rnoise_RUN_stim_std(ist).*[-1; 1],...
        '-','Color',colors(allStim(ist),:),'LineWidth',4)
    plot(allStim(ist),Rnoise_RUN_stim(ist),'.','Color',colors(allStim(ist),:),'MarkerSize',50)
    plot(allStim(ist), Rnoise_stim(ist) ,'x','Color',colors(allStim(ist),:),'MarkerSize',20,'LineWidth',4)
end

set(gca,'xtick',allStim,'xticklabel',Info.stim_ID_key(allStim)',...
    'xlim',[0 ist+1])

title(sprintf('Clus: %i (%s) & %i (%s)',Clu1,UnitInfo.SpkShape{Un1},Clu2,UnitInfo.SpkShape{Un2} ))

print_eps_kp(gcf,fullfile(savedir,sprintf('AllStim_%s_%s_Clus_%i_%i_500msConv',SUBJECT,SESSION,Clu1,Clu2)))



end