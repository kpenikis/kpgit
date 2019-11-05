function allRastersSession_VS(SUBJECT,SESSION)
%
%  allRastersSession(SUBJECT, SESSION )
%    Plots a raster and psth for each stimulus, for all SU from the session.
%    Uses Clusters struct, not UnitData table. 
%
%  KP, 2019-01, updated 2019-02
%


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

q = load(fullfile(fn.processed,'UnitsVS'));
UnitData  = q.UnitData;
UnitInfo  = q.UnitInfo;
clear q

filename  = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename  = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename  = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));


% Filter Unit files to just this session and sort by baseline FR
UnitData  = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo  = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
[~, baseFRrank] = sort([UnitData.BaseFR]);
UnitInfo  = UnitInfo(baseFRrank,:);
UnitData  = UnitData(baseFRrank);
Clusters  = Clusters(baseFRrank);

% Get indices of narrow spikes and regular spiking units
UnitInfo = labelRSNS(UnitInfo);


%% Prepare figures

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
rng('shuffle')

scrsz = get(0,'ScreenSize');  %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/6 scrsz(4)];

figwidthscales = UnitData(1).Dur./1000;


% Set colors
colors = [ ...
      0 200 150;...
      0 200 150;...
      0   0   0;...
      0   0   0;...
    181   0  52;...
    255  87  51]./255;

rasterdotsize  = 5;
psthlinewidth  = 4;

raster_colors = [0.2 0.2 0.2; 0.5 0.5 0.5];


%%

Stimuli = unique([TrialData.trID])';
Stimuli(Stimuli==0)=[];

AvgEnv = nan(numel(Stimuli),6000);
PSTH   = nan(numel(Stimuli),6000);

for ist = Stimuli
    
    stid = Stimuli(ist);
    
    add_y = 0;%zeros(1,7);
    N_un  = 0;%zeros(1,7);
    
    %%
    % Set up figure
    hf(ist) = figure;
    set(gcf,'Position',tallrect.* [1 1 figwidthscales(stid) 1])
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
    hs(ist,2) = subplot(9,1,2:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'ytick',[],...
        'xticklabel',[0 (figwidthscales(stid)*1000)])
    
    theseClus = 1:numel(Clusters); %[13 4 11 12 13 15]
    for iUn = theseClus 
        
        % Get spiketimes
        thisClu = Clusters(iUn);
        maxChan = thisClu.maxChannel;
        %originally: seconds, not rounded
        spiketimes = unique(round(thisClu.spikeTimes*1000)');
        
        % Get stimulus params
        dBSPL = UnitData(iUn).spl;
        LP    = UnitData(iUn).lpn;
        
        % Convert FR to z-score
        bs_smth = 20;
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        
        
        % Get all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(maxChan).trials,dBSPL,LP,0);
        
        
        %% Collect trial indices and timestamps
        
        TDidx = [];
        TDidx = all_TDidx([TrialData.trID(all_TDidx)]==stid);
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        
        t3 = t2 + Duration;
        
        
        %%
        % Preallocate
        raster_x = [];
        raster_y = [];
        stim   = nan( numel(TDidx), Duration+1 );
        psth   = nan( numel(TDidx), Duration+1 );
        
        % Collect spikes/FR/rms for this stimulus/unit
        for it = 1:numel(TDidx)
            
            stim(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) - 1;
            
            raster_x = [raster_x sp];
            raster_y = [raster_y it*ones(1,numel(sp))];
            
            psth(it,:) = Stream_FRsmooth(1, t2(it) : t3(it) );
            
        end %it
        
        
        %% Add to plots
        
        figure(hf(ist)); hold on
        
        % Stimulus RMS
        subplot(hs(ist,1)); hold on
        plot(0:Duration, mean(stim,1,'omitnan'),...
            'LineWidth',psthlinewidth,'Color',colors(stid,:))
        set(gca,'Color','none','xtick',[],'ytick',[])
        
        % Raster
        subplot(hs(ist,2)); hold on
        plot(raster_x, raster_y + add_y(end),...
            '.','MarkerSize',rasterdotsize,'Color',raster_colors(1+mod(N_un,2),:))
        set(gca,'Color','none','xtick',[],...
            'tickdir','out','ticklength',40/Duration.*[1 1])
        box off
        
        
        add_y = [add_y add_y(end) + numel(TDidx)];
        N_un = N_un + 1;
        
        
        % Save data for shihab
        AvgEnv(ist,1:size(stim,2)) = mean(stim,1);
        PSTH(ist,1:size(psth,2))   = mean(psth,1);
%         savedataname = sprintf('VSresponse_Clu%i_Stim%i',thisClu.clusterID,ist);
%         save(fullfile(fn.root,'SharedData',savedataname),'AvgEnv','PSTH','-v7.3')
        
    end % iUn
    
    if add_y(end) == 0
        continue
    end
    
    suptitle(sprintf('"%s"\n%s %s\n%i units',...
        Info.stim_ID_key{stid}, SUBJECT,SESSION,N_un ))
    
    CluLabels = cellfun(@num2str,num2cell([UnitInfo(theseClus,:).Clu]),'UniformOutput',false)';
    set(gca,'ytick',add_y(2:end),'yticklabel',CluLabels,'ylim',[0 add_y(end)+1],...
        'xtick',0:500:Duration,'xticklabel',0:500:Duration)
    xlabel('Time (ms)')
    ylabel([num2str(mode(diff(add_y))) ' trials each'])
    
    subplot(hs(ist,1));
    set(gca,'ylim',[0 1])
    
    
    %% Save figure
    
    savedir = fullfile(fn.figs,'VocodedSpeech','SessionRasters',sprintf('%s_%s',SUBJECT,SESSION));
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    savename = sprintf('%s_%s_%s',SUBJECT,SESSION,Info.stim_ID_key{ist});
    print_eps_kp(hf(ist),fullfile(savedir,savename))
    
    
end %ist

keyboard


% Save data for shihab

AvgEnv = (AvgEnv-min(min(AvgEnv)));
AvgEnv = AvgEnv / max(max(AvgEnv));

Info.stim_ID_key{3} = 'ICantBlabSuch';
Info.stim_ID_key{4} = 'ImTheLorax';
StimInfo=Info;
Info=struct();
Info.fs             = 1000;
Info.rows_stim_label = StimInfo.stim_ID_key;

savename = sprintf('VSresponses_Clu%i', thisClu.clusterID);
save(fullfile(fn.root,'SharedData',savename),'AvgEnv','PSTH','Info','-v7.3')



end