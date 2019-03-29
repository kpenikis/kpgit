function plotRasters_allUnits(stimID)
%
%  plotRasters_allUnits()
%    Plots a raster and psth for each stimulus.
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2018-07
%



%% Load Unit files

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


N = 0; %zeros(1,7);
add_trs = 0; %zeros(1,7);
ymaxval = 50;
rng('shuffle')
ist = 1;


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',20)

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937];

% Set colors (stimuli)
% colors = [   0 200 150;...
%             84  24  69;...
%            120  10  41;...
%            181   0  52;...
%            255  87  51;...
%            255 153   0]./255;
% colors = [ colors; ...
%             [37  84 156]./255 ;...
%             [19 125 124]./255 ];

rasterdotsize  = 6;
psthlinewidth  = 0.75;


% Set colors (unit response types)
colors = [   0 200 150;...
    84  24  69;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;
colors = [ colors; ...
    [37  84 156]./255 ;...
    [19 125 124]./255 ];

raster_colors = [0.3 0.7 0.7; 0.8 0.8 0.3];
bkgrnd_colors = [0.8  0   0 ;  0   0  0.8];


% Set up figure
hf = figure;
set(gcf,'Position',halfscreen.* [1 1 figwidthscales(stimID) 1])
hold on
hs(1,1) = subplot(9,1,1);
set(gca,'xlim',[0 figwidthscales(stimID)*1000],'xtick',[],'ytick',[])
hs(1,2) = subplot(9,1,2:6);
set(gca,'xlim',[0 figwidthscales(stimID)*1000],'xtick',[],'ytick',[])
hs(1,3) = subplot(9,1,7:9);
set(gca,'xlim',[0 figwidthscales(stimID)*1000],'xtick',[0 (figwidthscales(stimID)*1000)],'xticklabel',[0 (figwidthscales(stimID)*1000)])



%%

% Sort units by overall avg FR resp to all stim
yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {UnitData.FR_raw_tr},'UniformOutput',false);
Units_meanFRs = mean(vertcat(yyy{:}),2,'omitnan');

[~,isrt] = sort(Units_meanFRs);
Units_meanFRs = Units_meanFRs(isrt);
Unit_sort = UnitData(isrt);


% Identify units labeled AM responsive
if ~isfield(UnitData,'iBMF_FR')
    keyboard
    [sigUnits,Unit_sort] = identifyResponsiveUnits(Unit_sort);
end


theseUnits = 1:numel(Unit_sort);


for iUn = theseUnits
    
    % Get this unit's info
    subject = Unit_sort(iUn).Subject;
    session = Unit_sort(iUn).Session;
    channel = Unit_sort(iUn).Channel;
    clu     = Unit_sort(iUn).Clu;
    
    if numel(session)>2, continue, end
    
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    if (iUn>1 && ~( strcmp(subject,Unit_sort(iUn-1).Subject) && strcmp(session,Unit_sort(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,Unit_sort(iUn-1).Subject) && strcmp(session,Unit_sort(iUn-1).Session) && channel==Unit_sort(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % Get stimulus params
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    % Get spiketimes
    spikes = Spikes.sorted(channel);
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    
    
    % Convert FR to z-score
    bs_smth = 20;
    Stream_FRsmooth = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    allStim = unique(TrialData.trID(all_TDidx))';
    
    if ~any(allStim==stimID)
        continue
    end
    
    minTrs = min(Ntrials(~isnan(Ntrials)));
    
    % Skip if too few trials
    if minTrs<20
        continue
    elseif minTrs>30
        minTrs = 20;
    end
    
    N = N+1;
    
    
    
    %% Collect trial indices and timestamps
    
    if stimID==3 || stimID==6
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stimID & TrialData.ITIflag(all_TDidx)==0 );
        % Find Pdc trials that follow same rate during ITI
        TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stimID);
        
        TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stimID & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
        TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        
    else
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stimID );
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
    
    
    
    kt = randperm(length(t2));
    kt = kt(1:minTrs);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = t2 + Duration;
    
    
    %%
    
    % Preallocate
    raster_x = [];
    raster_y = [];
    stim   = nan( numel(TDidx), Duration+1 );
    psth   = nan( numel(TDidx), Duration+1 );
    
    
    % Collect spikes/FR/rms for this stimulus/unit
    for it = 1:numel(TDidx)
        
        psth(it,:) = ...
            Stream_FRsmooth( t2(it) : t3(it) );
        
        stim(it,:) = ...
            SoundStream(1, t2(it) : t3(it) )...
            ./ max(SoundStream(1, t2(it) : t3(it) ));
        
        sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
            & spiketimes<=t3(it) ) - t2(it) - 1;
        
        raster_x = [raster_x sp];
        raster_y = [raster_y it*ones(1,numel(sp))];
        
    end %it
    
    
    
    %% Add to plots
    
    figure(hf(ist)); hold on
    
    % Stimulus
    subplot(hs(ist,1)); hold on
    plot(0:Duration, mean(stim,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','k')
    set(gca,'Color','none','xtick',[],'ytick',[])
    hold off
    
    % Raster
    
    subplot(hs(ist,2)); hold on
    
    
    plot(raster_x, raster_y + add_trs,...
        '.','MarkerSize',rasterdotsize,'Color',raster_colors(1+mod(N,size(colors,1)),:))
    
    add_trs = add_trs + max(raster_y);
    
    set(gca,'Color','none','xtick',[],'ytick',[],'ylim',[0 add_trs+1],...
        'tickdir','out','ticklength',40/Duration.*[1 1])
    box off
    hold off
    
    
    % PSTH
    subplot(hs(ist,3)); hold on
%     plot(0:Duration, mean(psth,1,'omitnan'),...
%         'LineWidth',psthlinewidth,'Color',colors(1+mod(N,size(colors,1)-1),:))
%     set(gca,'Color','none',...
%         'tickdir','out','ticklength',40/Duration.*[1 1])
    hold off
    
    ymaxval = max(ymaxval,max(mean(psth,1,'omitnan')));
    
    
    
end % iUn


% Finish figure

set(hs(:,3),'ylim',[0 ceil(ymaxval)],'ytick',ceil(ymaxval))
ylabel('Spikes/sec')
subplot(hs(ist,2)); hold on
ylabel('20-30 trs per unit')

suptitle(sprintf('%s  |  N = %i units with >=20 trials', Info.stim_ID_key{stimID}, N ))



%% Save figure

savedir = fullfile(fn.processed,'Rasters');
if ~exist(savedir,'dir')
    mkdir(savedir)
    mkdir([savedir '/eps'])
    mkdir([savedir '/svg'])
end

savename = sprintf('AllUnits_%s',Info.stim_ID_key{stimID});

print_eps_kp(hf(ist),fullfile([savedir '/eps'],savename))
print_svg_kp(hf(ist),fullfile([savedir '/svg'],savename))




end