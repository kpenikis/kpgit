function screenUnitsRasters(SUBJECT,SESSION)
%
%  screenRasters(SUBJECT, SESSION )
%    Plots a raster and psth for each stimulus, for all SU from the session.
%    Excludes datapoints based on: min Ntrials, min FR.
%
%  KP, 2019-01
%


UseTempFolder = 1;


%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));


% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];



%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
rng('shuffle')

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

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


%%
for iUn = 1:numel(Clusters)    
    
    close all
    
    % FOR NOW, ADD PLACEHOLDER FOR ARTIFACT TRIALS FOR SYNAPSE/KS DATA
    if ~isfield(Info,'artifact')
        Info.artifact(64).trials = [];
    end
    if ~isfield(Info,'stim_ID_key')
        Info.stim_ID_key = { 'Warn';  '2';   '4';   '8';  '16';  '32';   'AC';  'DB'  };
    end
    
    
    % Get spiketimes
    thisClu = Clusters(iUn); 
    maxChan = thisClu.maxChannel;
    %originally: seconds, not rounded
    spiketimes = round(thisClu.spikeTimes*1000)';
    
    
    % Get stimulus params
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    
    % Convert FR to z-score
    bs_smth = 20;
    [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    meanFR = 1000* sum(Stream_Spikes(TrialData.onset(1):end)) / length(Stream_Spikes(TrialData.onset(1):end));
    
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(maxChan).trials,dBSPL,LP);
    allStim = unique(TrialData.trID(all_TDidx))';
    
    minTrs = min(Ntrials(~isnan(Ntrials)));
    
    % Skip if too few trials
%     if minTrs<10, continue, end
        
    ymaxval = 0;

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
                
        kt = randperm(length(t2));
        kt = kt(1:minTrs);
%         kt = 1:length(t2);
        t2     = t2(kt);
        TDidx  = TDidx(kt);
        
        t3 = t2 + Duration;
                
        
        %% 
        % Set up figure
        hf(ist) = figure;
        set(gcf,'Position',largerect.* [1 1 figwidthscales(stid) 1])
        hold on
        hs(ist,1) = subplot(9,1,1);
        set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
        hs(ist,2) = subplot(9,1,2:6);
        set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[],...
            'ylim', [0 length(t2)+1])
        hs(ist,3) = subplot(9,1,7:9);
        set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'xticklabel',[0 (figwidthscales(stid)*1000)])
        
        
        
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
            'LineWidth',psthlinewidth,'Color',colors(stid,:))
        set(gca,'Color','none','xtick',[],'ytick',[])
        hold off
        
        % Raster
        if numel(raster_x)>0
            
        subplot(hs(ist,2)); hold on
            plot([raster_x; raster_x], [raster_y; raster_y] + 0.45.*[-ones(size(raster_y)); ones(size(raster_y))],...
                'Color',colors(stid,:),'LineWidth',3)
            set(gca,'Color','none','xtick',[],'ytick',it,'yticklabel',it,'ylim',[0 it+1],...
                'tickdir','out','ticklength',40/Duration.*[1 1])
            box off
            hold off
            
        end
        
        % PSTH
        subplot(hs(ist,3)); hold on
        plot(0:Duration, mean(psth,1,'omitnan'),...
            'LineWidth',psthlinewidth,'Color',colors(stid,:))
        set(gca,'Color','none',...
            'tickdir','out','ticklength',40/Duration.*[1 1])
        hold off
        
        ymaxval = max(ymaxval,max(mean(psth,1,'omitnan')));
        
        suptitle(sprintf('%s\n%s %s\nShank %i | maxCh: %i | cluID: %i',...
            Info.stim_ID_key{stid}, SUBJECT,SESSION,thisClu.shank,maxChan,thisClu.clusterID ))
        
        
    end %ist
    
    linkaxes(hs(:,3),'y')
    set(hs(:,3),'ylim',[0 ceil(ymaxval)],'ytick',ceil(ymaxval) )
    
    
    %% Plot unit summary
    
    hfClu = figure;
    set(gcf,'Position',largerect.* [1 1 1.5 0.5])
    hold on
    subplot(1,3,1);
    plot(thisClu.maxChTemp,'k','LineWidth',3)
    set(gca,'xlim',[1 length(thisClu.maxChTemp)],'xtick',[],'Color','none')
    title(sprintf('%i events | %0.1f Hz',length(spiketimes),meanFR ))
    subplot(1,3,2);
    histogram([-diff(spiketimes) diff(spiketimes)],[-60:60],'EdgeColor','none','FaceColor','k','FaceAlpha',1);
    set(gca,'xlim',[-50 50],'Color','none')
    title(sprintf('%0.1f%% violation rate', 100*sum(diff(spiketimes)<3)/length(spiketimes) ))
    subplot(1,3,3);
    try
    histogram([thisClu.Amplitudes],50,'FaceColor','k','EdgeColor','none')
    end
    title('Distribution of amplitudes')
    
    suptitle(sprintf('Shank %i  |  maxCh: %i  |  cluID: %i',thisClu.shank,maxChan,thisClu.clusterID ))
    
    
    %% Save figures
    
    if UseTempFolder
%         savedir = fullfile(fn.processed,'Rasters_TEMP',sprintf('%s_%s_%i_%i',SUBJECT,SESSION,thisClu.shank,thisClu.clusterID));
        savedir = fullfile(fn.sessdata,'Rasters','tmp',sprintf('%s_%s_%i_%i',SUBJECT,SESSION,thisClu.shank,thisClu.clusterID));
    else
        savedir = fullfile(fn.sessdata,'Rasters',sprintf('%s_%s_%i_%i',SUBJECT,SESSION,thisClu.shank,thisClu.clusterID));
    end
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    for ist = 1:numel(allStim)
        savename = sprintf('%s_%s_%i_%i_clu%i_%s',SUBJECT,SESSION,maxChan,thisClu.shank,thisClu.clusterID,Info.stim_ID_key{ist});
        print_eps_kp(hf(ist),fullfile(savedir,savename))
    end
    
    savename = sprintf('%s_%s_%i_%i_clu%i',SUBJECT,SESSION,thisClu.shank,maxChan,thisClu.clusterID);
    print_eps_kp(hfClu,fullfile(savedir,savename))
    
    close(hf); close(hfClu)
    clear hf hs hfClu
    
end % iUn

% keyboard





end