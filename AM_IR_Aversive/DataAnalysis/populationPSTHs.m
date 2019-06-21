function populationPSTHs(SUBJECT,SESSION)
% 
% Must be single session.
% Here, plot aggregate rasters and histograms, combining narrow spiking and
% regular spiking units.
% 
% KP, 2019-03-15
% 


% close all
global fn 


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

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

% Filter Unit files to just this session and sort by baseline FR
UnitData = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
[~, baseFRrank] = sort([UnitData.BaseFR]);
UnitInfo = UnitInfo(baseFRrank,:);
UnitData = UnitData(baseFRrank);

% Get indices of narrow spikes and regular spiking units
iNarrow = find([UnitInfo.WidthHalfMax]<=0.23 & [UnitInfo.TroughPeak]<0.5)';
iRegSpk = find([UnitInfo.WidthHalfMax]>0.23 & [UnitInfo.TroughPeak]>=0.5)';


% Also find trials to skip 
Channels = unique([UnitInfo.Channel]);
artifactTrs = [];
for ich = Channels'
    artifactTrs = [artifactTrs Info.artifact(ich).trials'];
end
artifactTrs = unique(artifactTrs);


%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallfig  = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];

% Set colors
ColNarrow = [0.92 0.24 0.37];
ColRegSpk = [0.11 0.16 0.55];

        
% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];

minTrs = 12;


%% 

for stid = unique([TrialData.trID])'
    
    if stid==0, continue, end
    
    % Prepare the figure
    hf(stid) = figure;
    set(hf(stid),'Position',tallfig.* [1 1 figwidthscales(stid) 1],'NextPlot','add')
    hold on
    hs(1)=subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[],'Color','none');
    box off
    hs(2)=subplot(9,1,2:3);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'Color','none');
    box off
    hs(3)=subplot(9,1,4:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'ytick',[],...
        'xticklabel',[0 (figwidthscales(stid)*1000)],'Color','none','tickdir','out');
    xlabel('Time (ms)')
    box off
    
    
    % Preallocate
    N_un = 0;
    PSTH_Nar = [];
    PSTH_Reg = [];
    
    for iUn = 1:size(UnitInfo,1)
        
        if ismember(iUn,iNarrow)
            unCol = ColNarrow;
        elseif ismember(iUn,iRegSpk)
            unCol = ColRegSpk;
        else
            continue
        end
        
        channel = UnitData(iUn).Channel(1);
        clu     = UnitData(iUn).Clu(1);
        
        % Get spiketimes (KS)
        spiketimes = round(Clusters(([Clusters.maxChannel]==channel & [Clusters.clusterID]==clu)).spikeTimes*1000 - spkshift)';
        
        % Get sound parameters
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,50,'silence');
        
        
        %%
        % Find all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,artifactTrs,dBSPL,LP);
        
        allStim = unique(TrialData.trID(all_TDidx));
        
        if sum(Ntrials < minTrs)==1
            keyboard
            all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
            allStim(Ntrials<minTrs)  = [];
            Ntrials(Ntrials<minTrs) = [];
        elseif  sum(Ntrials < minTrs)>1
            keyboard
        end
        
        %%
        % Get trial numbers for this stimulus
        
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
        
        
        % Timestamps of onsets and offsets
        
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        
        % Add ITI trials (shortened to match duration)
        if ~isempty(TDidx_iti)
            t2 = [t2; TrialData.onset(TDidx_iti)];
            TDidx = [TDidx; TDidx_iti];
        end
        
        t3 = t2 + Duration;
        
        
        % Collect spikes/FR/rms for this stimulus/unit
        raster_x = [];
        raster_y = [];
        stim   = nan( numel(TDidx), Duration+1 );
        psth   = nan( numel(TDidx), Duration+1 );
        
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
            
        end %it
        
        
        %% Add to plot
        
        subplot(hs(3)); hold on
        plot([raster_x; raster_x], [raster_y; raster_y] + 0.45.*[-ones(size(raster_y)); ones(size(raster_y))],...
            'Color',unCol,'LineWidth',2)
        
        N_un = N_un + 1;
        
        
        %% Aggregate PSTH
        
        if ismember(iUn,iNarrow)
            PSTH_Nar = [PSTH_Nar; mean(psth,1)];
        elseif ismember(iUn,iRegSpk)
            PSTH_Reg = [PSTH_Reg; mean(psth,1)];
        end
        
        
    end %iUn
    
    set(gca,'ylim',[0 it+0.5],'xlim',[0 Duration])
    ylabel([num2str(it) ' trials'])
    
    
    % Plot PSTHs
    subplot(hs(2)); hold on
    plot(sum(PSTH_Nar,1),'Color',ColNarrow,'LineWidth',3) %/max(mean(PSTH_Nar,1))
    plot(sum(PSTH_Reg,1),'Color',ColRegSpk,'LineWidth',3) %/max(mean(PSTH_Reg,1))
    set(gca,'ylim',[0 max([sum(PSTH_Nar,1) sum(PSTH_Reg,1)])])
    ylabel('Firing rate (sp/s)')
    
    % Plot stimulus RMS
    subplot(hs(1)); hold on
    plot(mean(stim,1),'k','LineWidth',3)
    
    suptitle(Info.stim_ID_key{stid})
    
    savedir = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/SessionPopulationRasters';
    savename = sprintf('%s_%s_%s',SUBJECT,SESSION,Info.stim_ID_key{stid});
    print_eps_kp(hf(stid),fullfile(savedir,savename))
    
    
    %% Normalized PSTHs
    
%     PSTH_Nar_norm = PSTH_Nar - repmat(min(PSTH_Nar,[],2),1,size(PSTH_Nar,2));
%     PSTH_Nar_norm = PSTH_Nar_norm ./ repmat(max(PSTH_Nar_norm,[],2),1,size(PSTH_Nar_norm,2));
%     
%     PSTH_Reg_norm = PSTH_Reg - repmat(min(PSTH_Reg,[],2),1,size(PSTH_Reg,2));
%     PSTH_Reg_norm = PSTH_Reg_norm ./ repmat(max(PSTH_Reg_norm,[],2),1,size(PSTH_Reg_norm,2));
%     
%     hf2(stid) = figure;
%     set(hf2(stid),'Position',tallfig.* [1 1 figwidthscales(stid) 1],'NextPlot','add')
%     hold on
%     imagesc([PSTH_Reg_norm; PSTH_Nar_norm])
%     colormap(cmocean('ice'))
%     xlim([0 Duration])
%     ylim(0.5+[0 N_un])
%     
%     title(Info.stim_ID_key{stid})
    
    
end %stid




end