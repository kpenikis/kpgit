function plotRasters_eachUnit(opt)
%
%  plotRasters_eachUnit(subject, session, [channel, clu] )
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

N=0;
rng('shuffle')


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',20)

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937];

add_y = zeros(1,7);
N_un  = zeros(1,7);

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

% Sort units by overall avg FR resp to all stim
yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {UnitData.FR_raw_tr},'UniformOutput',false);
Units_meanFRs = mean(vertcat(yyy{:}),2,'omitnan');

[~,isrt] = sort(Units_meanFRs);
Units_meanFRs = Units_meanFRs(isrt);
Unit_sort = UnitData(isrt);


% Identify units labeled AM responsive
if ~isfield(UnitData,'iBMF_FR')
%     keyboard
    [sigUnits,Unit_sort] = identifyResponsiveUnits(Unit_sort);
end

if nargin>0
    theseUnits = find(strcmp({Unit_sort.Session},'Sep17-AM'));
%         theseUnits = theseUnits(1);
else
    theseUnits = 1:numel(Unit_sort);
end

for iUn = theseUnits
    
    % Get this unit's info
    subject = Unit_sort(iUn).Subject;
    session = Unit_sort(iUn).Session;
    channel = Unit_sort(iUn).Channel;
    clu     = Unit_sort(iUn).Clu;
        
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    
%     loadThisUnitData(subject,session)
    if (iUn>1 && ~( strcmp(subject,Unit_sort(iUn-1).Subject) && strcmp(session,Unit_sort(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        clear Info TrialData
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,Unit_sort(iUn-1).Subject) && strcmp(session,Unit_sort(iUn-1).Session) && channel==Unit_sort(iUn-1).Channel ) )  || iUn==1
        clear Spikes Clusters 
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % FOR NOW, ADD PLACEHOLDER FOR ARTIFACT TRIALS FOR SYNAPSE/KS DATA
    if ~isfield(Info,'artifact')
        Info.artifact(64).trials = [];
    end
    Info.stim_ID_key = { 'Warn';  '2';   '4';   '8';  '16';  '32';   'AC';  'DB'  };
    
    
    % Get spiketimes
    if exist('Spikes','var')
        spikes = Spikes.sorted(channel);
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    elseif exist('Clusters','var')
        spiketimes = round(Clusters([Clusters.clusterID]==clu & [Clusters.maxChannel]==channel).spikeTimes *1000)';
    end
    
    
    % Get stimulus params
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    
    % Convert FR to z-score
    bs_smth = 20;
    [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    allStim = unique(TrialData.trID(all_TDidx))';
    
    minTrs = min(Ntrials(~isnan(Ntrials)));
    
    % Skip if too few trials
%     if minTrs<10, continue, end
    
    STIMidx = find(allStim==2 | allStim==7 | allStim==8);
    
    close all
    ymaxval = 0;

    for ist = 1:numel(allStim)
        
        stid = allStim(ist);
        
        
%         if nargin<1
            % Set up figure
            hf(ist) = figure;
            set(gcf,'Position',largerect.* [1 1 figwidthscales(stid) 1])
            hold on
            hs(ist,1) = subplot(9,1,1);
            set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
            hs(ist,2) = subplot(9,1,2:6);
            set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
            hs(ist,3) = subplot(9,1,7:9);
            set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'xticklabel',[0 (figwidthscales(stid)*1000)])
%         end
        
        
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
            'LineWidth',psthlinewidth,'Color',colors(stid,:))
        set(gca,'Color','none','xtick',[],'ytick',[])
        hold off
        
        % Raster
        if numel(raster_x)>0
            
        subplot(hs(ist,2)); hold on
            plot(raster_x, raster_y,...
                '.','MarkerSize',rasterdotsize,'Color',colors(stid,:))
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
        
        try
            suptitle(sprintf('%s   |   %s %s ch%i clu%i\nN=%i, avg FR = %0.2f Hz  |  BMF = %i Hz ',...
                Info.stim_ID_key{stid}, subject,session,channel,clu,iUn,Units_meanFRs(iUn),AMrates(Unit_sort(iUn).iBMF_FR) ))
        catch
            suptitle(sprintf('%s   |   %s %s ch%i clu%i\nN=%i, avg FR = %0.2f Hz',...
            Info.stim_ID_key{stid}, subject,session,channel,clu,iUn,Units_meanFRs(iUn) ))
        end        
        
        
    end %ist
    
    
    %% Save figure
    
    linkaxes(hs(:,3),'y')
    set(hs(:,3),'ylim',[0 ceil(ymaxval)],'ytick',ceil(ymaxval) )
    
    savedir = fullfile(fn.processed,'Rasters',sprintf('%s_%s_%i_%i',subject,session,channel,clu));
    if ~exist(savedir,'dir')
        mkdir(savedir)
        mkdir([savedir '/eps'])
%         mkdir([savedir '/svg'])
    end
    
    for ist = 1:numel(allStim)
        
        savename = sprintf('%s_%s_%i_%i_%s',subject,session,channel,clu,Info.stim_ID_key{ist});
        
        print_eps_kp(hf(ist),fullfile([savedir '/eps'],savename))
%         print_svg_kp(hf(ist),fullfile([savedir '/svg'],savename))
    end
    
    close(hf)
    clear hf hs
    
end % iUn

% keyboard





end