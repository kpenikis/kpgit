function [FR_allSt,std_allSt] = PlotExampleResponses(SUBJECT,SESSION, CluID, maxChan)
%
%  PlotExampleResponses( SUBJECT, SESSION, CluID )
%    Plots a raster and psth for each stimulus for the designated cell.
%
%   NEXT: add tuning curves for VS and trial variability
%
%  KP, 2019-01, updated 2019-02
%

global AMrates minTrs

%%%%%%%%%%%%%%%%%%%%
TTP = 20;
%%%%%%%%%%%%%%%%%%%%
minTrs = 10;
%%%%%%%%%%%%%%%%%%%%
% USE_MEASURE = 'FR';
%%%%%%%%%%%%%%%%%%%%
AMrates = [2 4 8 16 32];
%%%%%%%%%%%%%%%%%%%%

%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

iUn = find(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT) & UnitInfo.Clu==CluID);
if numel(iUn)~=1, keyboard, end

% Get spiketimes
if nargin<4
    maxChan = UnitData(iUn).Channel;
    iClu = find([Clusters.clusterID] == CluID);
    spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
else
    spiketimes = unique(round( Spikes.sorted(maxChan).spiketimes(Spikes.sorted(maxChan).assigns==CluID') * 1000 - spkshift ));  %ms
end


%% Prepare figures

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',32)
psthlinewidth  = 4;
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
tcsquare   = [1 scrsz(4) scrsz(3)/3 scrsz(4)/2];
widescreen = [1 scrsz(4) scrsz(3) scrsz(4)/2];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];
Info.stim_ID_key{9} = 'Silence';
ymaxval = 0;

savedir = fullfile(fn.figs,'Rasters',sprintf('Un%i_%s_%s_%i',iUn,SUBJECT,SESSION,CluID));
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% Set colors
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];

%%

Stimuli   = [9 1:8];
FRtrials  = cell(1,numel(Stimuli));
FR_allSt  = nan(1,numel(Stimuli));
std_allSt = nan(1,numel(Stimuli));
FF_allSt  = nan(1,numel(Stimuli));

AvgEnv = nan(numel(Stimuli),2000);
PSTH   = nan(numel(Stimuli),2000);

for ist = Stimuli
    
    stid = ist;
    
    %%
    % Set up figure
    hf(ist) = figure;
    set(gcf,'Position',tallrect.* [1 1 figwidthscales(stid) 1])
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ylim',[0 1],'ytick',[])
    hs(ist,2) = subplot(9,1,2:5);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[])
    hs(ist,3) = subplot(9,1,6:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],...
        'xticklabel',[0 (figwidthscales(stid)*1000)])
    
    % Get stimulus params
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    % Convert FR to z-score
    bs_smth = 20;
    [Stream_FRsmooth,Stream_zscore,Stream_Spikes,~] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(maxChan).trials,dBSPL,LP,1);
    allStim = unique(TrialData.trID(all_TDidx))';
    
    Silence = false;
    if stid==9
        Silence = true;
    elseif ~ismember(stid,allStim)
        continue
    end
    
    
    %% Collect trial indices and timestamps
    
    if ~Silence                                          % SOUND TRIALS
        
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
        
    else                                               % SILENCE TRIALS
        clear t2 t3 Duration t_win
        SilPd = [TrialData.onset(1) TrialData.offset(1)];
        Duration = 1000;
        
        t2 = TrialData.onset(1) : Duration : (TrialData.offset(1)-mod(diff(SilPd),Duration)-Duration);
        TDidx = 1:length(t2);
        Ntrials(ist) = length(t2);
        
        TTP = min(min(Ntrials(Ntrials>0),40));
    end
    
%     kt     = randperm(Ntrials(ist),min(Ntrials));
    kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    if Silence && t3(end)>TrialData.offset(1)
        keyboard
    end
    
    
    %%
    % Preallocate
    raster_x = [];
    raster_y = [];
    stim   = nan( numel(TDidx), Duration+1 );
    psth   = nan( numel(TDidx), Duration+1 );
    FR     = nan( numel(TDidx), 1 );
    
    % Collect spikes/FR/rms for this stimulus/unit
    for it = 1:numel(TDidx)
        
        stim(it,:) = ...
            SoundStream(1, t2(it) : t3(it) )...
            ./ max(SoundStream(1, t2(it) : t3(it) ));
        
        sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
            & spiketimes<=t3(it) ) - t2(it) - 1;
        
        raster_x = [raster_x sp];
        raster_y = [raster_y it*ones(1,numel(sp))];
        
        psth(it,:) = ...
            Stream_FRsmooth( t2(it) : t3(it) );
        
        FR(it)   = numel(sp)/(Duration/1000);
        
    end %it
    
    if Silence
        stim = zeros(size(stim));
    end
    
    %% Save some data
    
    FRtrials{ist} = FR;
    
    % FR and std
    FR_allSt(ist)  = mean(FR);
    std_allSt(ist) = std(FR);
    
    % FF
    FF_allSt(ist)  = var(FR)/mean(FR);
    
    
    %% Add to plots
    
    figure(hf(ist)); hold on
    
    %--------- RMS ---------
    subplot(hs(ist,1)); hold on
    plot(0:Duration, mean(stim,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','k')
    box off
    set(gca,'Color','none','xtick',[],'ytick',[],'visible','off')    
    title(sprintf('%s    .    %s    .    %s    .    Clu %i    .    iUn=%i',...
        Info.stim_ID_key{stid}, SUBJECT,SESSION,CluID,iUn ),'FontSize',18)
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    
    %-------- Raster --------
    subplot(hs(ist,2)); hold on
    if ~isempty(raster_y(raster_y<=TTP))
        plot([raster_x(raster_y<=TTP); raster_x(raster_y<=TTP)], raster_y(raster_y<=TTP) + [-0.5; 0.5] ,...
            'Color','k','LineWidth',psthlinewidth/2)
    end
    set(gca,'Color','none','xtick',[])
    box off
    
    set(gca,'ytick',[],'ylim',[0.5 TTP+0.5],'xtick',[],'visible','off')
    ylabel(sprintf('%i of %i trials',TTP,length(t2)))
    
    %--------- PSTH ---------
    subplot(hs(ist,3)); hold on
    plot(0:Duration, mean(psth,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','k')
    set(gca,'Color','none','xtick',unique([0:500:Duration Duration]),...
        'xticklabel',unique([0:500:Duration Duration]),'tickdir','out','ticklength',0.01.*[1 1])
    xlabel('Time (ms)')
    box off
    hold off
    
    ymaxval = max(ymaxval,max(mean(psth,1,'omitnan')));
    
    % Save RMS and PSTH data
    AvgEnv(ist,1:size(stim,2)) = mean(stim,1,'omitnan');
    PSTH(ist,1:size(psth,2))   = mean(psth,1,'omitnan');

    
end %ist


% Save data for shihab

% AvgEnv = (AvgEnv-min(min(AvgEnv)));
% AvgEnv = AvgEnv / max(max(AvgEnv));
% AvgEnv = AvgEnv(2:end,:);
% 
% PSTH   = PSTH(2:end,:);
% 
% StimInfo=Info;
% Info=struct();
% Info.fs             = 1000;
% Info.rows_stim_label = StimInfo.stim_ID_key(2:end-1);
% 
% savename = sprintf('AMresponses_Clu%i', CluID);
% save(fullfile(savedir,savename),'AvgEnv','PSTH','Info','-v7.3')


[IR_Prediction_sim, IR_Pred_std_sim, pvals] = predictIRresponse_simulation(FRtrials);



%% Tuning curves 

htc=figure;
set(gcf,'Position',widescreen)
hold on

xvals = [1:5 allStim(allStim>6)+2]; 


% - - - - - - - - - - - - - -   FR   - - - - - - - - - - - - - -

y_vals = FR_allSt;
y_errs = std_allSt;

subplot(1,3,1);
hold on

% Plot baseline FR
if abs(FR_allSt(9) - UnitData(iUn).BaseFR) > (FR_allSt(9)/20)
    keyboard
end
plot([0 max(xvals)+1],[FR_allSt(9) FR_allSt(9)],'--k','LineWidth',3)

% Plot periodic stimuli
plot(1:5,y_vals(2:6),'k','LineWidth',2)
for ir = 1:5
    plot([xvals(ir) xvals(ir)], y_vals(ir+1) + y_errs(ir+1)*[-1 1], ...
        '-','Color','k','LineWidth',4)
    plot(xvals(ir),y_vals(ir+1), 'o','MarkerSize',20,...
        'MarkerFaceColor','k','MarkerEdgeColor','none')
end

% Plot IR prediction
plot([7 7], IR_Prediction_sim + IR_Pred_std_sim*[-1 1],'k','LineWidth',4)
plot( 7, IR_Prediction_sim, 'ok','MarkerSize',20,...
    'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',4)

% Plot IR observed
for ir = 1:sum(xvals>5)
    plot(xvals(5+ir),y_vals(xvals(5+ir)-2),'o','MarkerSize',20,...
        'MarkerFaceColor',colors(xvals(5+ir)-2,:),'MarkerEdgeColor','none')
    plot([xvals(5+ir) xvals(5+ir)],y_vals(xvals(5+ir)-2) + y_errs(xvals(5+ir)-2)*[-1 1],...
        '-','Color',colors(xvals(5+ir)-2,:),'LineWidth',4)
    if pvals(xvals(5+ir)-8)<0.05
        plot(xvals(5+ir),y_vals(xvals(5+ir)-2) + 1.5*y_errs(xvals(5+ir)-2),'*k')
    end
end

% Finish formatting
set(gca,'XTick',sort([xvals 7]),...
    'XTickLabel',[Info.stim_ID_key(2:6)' 'Linear Prediction' Info.stim_ID_key(allStim(7:end))' ],...
    'TickLabelInterpreter','none')
xtickangle(45)

xlim([-1+min(xvals) max(xvals)+1])
keyboard
ylim([0 16])
ylabel('FR')
keyboard


% - - - - - - - - - - - - - -   VS   - - - - - - - - - - - - - -

subplot(1,3,2);
hold on

% Plot periodic stimuli
plot(1:5,UnitData(iUn).VSdata_spk(1,2:6),'-ok','MarkerSize',20,'LineWidth',2)
plot(find((UnitData(iUn).VSdata_spk(3,2:6)<0.001)),UnitData(iUn).VSdata_spk(1,1+find(UnitData(iUn).VSdata_spk(3,2:6)<0.001)),'.k','MarkerSize',55,'LineWidth',2)

% Finish formatting
set(gca,'XTick',sort([xvals 7]),...
    'XTickLabel',[Info.stim_ID_key(2:6)' 'Linear Prediction' Info.stim_ID_key(allStim(7:end))' ],...
    'TickLabelInterpreter','none')
xtickangle(45)

xlim([-1+min(xvals) max(xvals)+1])
ylim([0 1])
ylabel('VS')



% - - - - - - - - - - - - - -   FF   - - - - - - - - - - - - - -

subplot(1,3,3);
hold on

plot([0 max(xvals)+1],[FF_allSt(9) FF_allSt(9)],'--k','LineWidth',3)

% Plot periodic stimuli
plot(1:5,FF_allSt(2:6),'k','LineWidth',2)
for ir = 1:5
    plot(xvals(ir),FF_allSt(ir+1), 'o','MarkerSize',20,...
        'MarkerFaceColor','k','MarkerEdgeColor','none')
end

% Plot IR observed
for ir = 1:sum(xvals>5)
    plot(xvals(5+ir),FF_allSt(xvals(5+ir)-2),'o','MarkerSize',20,...
        'MarkerFaceColor',colors(xvals(5+ir)-2,:),'MarkerEdgeColor','none')
end

% Finish formatting
set(gca,'XTick',sort([xvals 7]),...
    'XTickLabel',[Info.stim_ID_key(2:6)' 'Linear Prediction' Info.stim_ID_key(allStim(7:end))' ],...
    'TickLabelInterpreter','none')
xtickangle(45)

xlim([-1+min(xvals) max(xvals)+1])
% keyboard
ylim([0 10])
ylabel('FF')
keyboard


%% Save figures

for ist = Stimuli
    
    figure(hf(ist));
    subplot(hs(ist,3));
    set(gca,'ylim',[0 ceil(ymaxval/10)*10],'ytick',linspace(0,ceil(ymaxval/10)*10,5))
    
    savename = sprintf('%s_%s_%i',Info.stim_ID_key{ist},SESSION,CluID);
    print_eps_kp(hf(ist),fullfile(savedir,savename))
    
end


% Save MTF figure
savename = sprintf('TuningCurves_%s_%i',SESSION,CluID);
print_eps_kp(htc,fullfile(savedir,savename))



end