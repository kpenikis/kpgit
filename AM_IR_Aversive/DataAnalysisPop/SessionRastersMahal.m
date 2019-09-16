function SessionRastersMahal(SUBJECT,SESSION)
%
%  SessionRastersMahal(SUBJECT, SESSION )
%    Plots a raster and psth for each stimulus, for all SU from the session.
%    Uses Clusters struct, not UnitData table. 
%    Based on allRastersSession.
% 
%  KP, 2019-08
%

% keyboard
% add TTP limit


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

% Sort units by baseline FR
[~, inspks] = sort(cellfun(@(x) sum(x>(TrialData.onset(1)/1000) & x<(TrialData.offset(1)/1000)), {Clusters.spikeTimes},'UniformOutput',true)); % baseline FR (during silence)
% [~, inspks] = sort(cellfun(@length, {Clusters.spikeTimes})); %overall nspikes
Clusters = Clusters(inspks);


%% Get mahal dist vector

gw_smth = 10;
Mdist = getPopMahalDistVector(Clusters,TrialData,spkshift,gw_smth);
q95  = quantile(Mdist,0.95);
q100 = max(Mdist);


%% Prepare figures

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
rasterdotsize  = 10;
psthlinewidth  = 4;
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];


% Set colors
colors = [   0 200 150;...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;
colors = [ colors; ...
    [37  84 156]./255 ;...
    [19 125 124]./255 ;...
    [ 0   0   0]./255];

raster_colors = [0 0 0; 0.2 0.5 0.8];
patch_colors  = [0.9 0.9 0.9; 1 1 1];

Info.stim_ID_key{9} = 'Silence';
rmsrange = [];


%%
Stimuli = 3:9;
for ist = Stimuli
    
    stid = Stimuli(ist);
    
    add_y    = 0;%zeros(1,7);
    N_un     = 0;%zeros(1,7);
    ytickstr = cell(1,numel(Clusters));
    
    %%
    % Set up figure
    hf(ist) = figure;
    set(gcf,'Position',tallrect.* [1 1 figwidthscales(stid) 1])
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
    hs(ist,2) = subplot(9,1,2:4);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
    hs(ist,3) = subplot(9,1,5:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'ytick',[],...
        'xticklabel',[0 (figwidthscales(stid)*1000)])
    
    
    Pop_psth   = nan( 300, figwidthscales(stid)*1000+1 , numel(Clusters));
    
    for iUn = 1:numel(Clusters)
        
        % Get spiketimes
        thisClu = Clusters(iUn);
        maxChan = thisClu.maxChannel;
        %originally: seconds, not rounded
        spiketimes = round(thisClu.spikeTimes*1000 - spkshift)';
        
        
        % Get stimulus params
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        % Convert FR to z-score
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,gw_smth,'silence');
        
        meanFR = 1000* sum(Stream_Spikes(TrialData.onset(1):end)) / length(Stream_Spikes(TrialData.onset(1):end));
        
        
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
            
            t2 = TrialData.onset(1) : Duration : (TrialData.offset(1)-mod(diff(SilPd),1000)-1000);
            TDidx = 1:length(t2);
            
        end
        
        
        kt = 1:length(t2);
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
        mdst   = nan( numel(TDidx), Duration+1 );
        
        % Collect spikes/FR/rms for this stimulus/unit
        for it = 1:numel(TDidx)
            
            stim(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
            psth(it,:) = ...
                Stream_FRsmooth(1, t2(it) : t3(it) );
            
            mdst(it,:) = ...
                Mdist(1, t2(it) : t3(it) )./max(Mdist).*ymaxval;
            
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) - 1;
            
            raster_x = [raster_x sp];
            raster_y = [raster_y it*ones(1,numel(sp))];
            
        end %it
        
        if Silence
            stim = zeros(size(stim));
        end
        
        Pop_psth(1:size(psth,1),:,iUn) = psth;
        
        
        %% Add to plots
        
        figure(hf(ist)); hold on
        
%         % Stimulus RMS
%         if iUn == numel(Clusters)
%             subplot(hs(ist,1)); hold on
%             plot(0:Duration, mean(stim,1,'omitnan'),...
%                 'LineWidth',psthlinewidth,'Color','k')
%             set(gca,'Color','none','xtick',[],'ytick',[])
%         end
        
        % Raster
        subplot(hs(ist,3)); hold on
        plot(raster_x, raster_y + add_y(end),...
            '.','Color',raster_colors(1+mod(N_un,2),:),'MarkerSize',rasterdotsize)
        set(gca,'Color','none','xtick',[],...
            'tickdir','out','ticklength',0.01.*[1 1])
        box off
        
        
        add_y = [add_y add_y(end) + max(kt)];
        N_un = N_un + 1;
        ytickstr{N_un} = num2str(thisClu.clusterID);
        
    end % iUn
    
    if add_y == 0
        continue
    end
    
    
    %% Collect data from population (stim, psth, and mahal)
    
    [Stream_FRsmooth,~,~,ymaxval] = convertSpiketimesToFR(round(vertcat(Clusters.spikeTimes)*1000-spkshift)',...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,gw_smth,'silence');
    
    % Collect data from trials (currently using last channel artifact info)
    stim   = nan( numel(TDidx), Duration+1 );
    psth   = nan( numel(TDidx), Duration+1 );
    mdst   = nan( numel(TDidx), Duration+1 );
    
    for it = 1:numel(TDidx)
        
        stim(it,:) = ...
            SoundStream(1, t2(it) : t3(it) )...
            ./ max(SoundStream(1, t2(it) : t3(it) ));
        
        psth(it,:) = ...
            Stream_FRsmooth(1, t2(it) : t3(it) );
        
        mdst(it,:) = ...
            Mdist(1, t2(it) : t3(it) )./q95.*ymaxval;
        
    end %it
    
    
    %% Quick check to compare population activity MEAN vs VARIANCE
    % cant inspect FF at single trial level without correcting for mean=0
    % no differences 
%     
%     ntr = find(sum(sum(~isnan(Pop_psth),2),3) == size(Pop_psth,2)*size(Pop_psth,3),1,'last');
%     
%     Pop_mean_tr = mean(Pop_psth(1:ntr,:,:),3);
%     Pop_var_tr  = var(Pop_psth(1:ntr,:,:),[],3);
%     Pop_ff_tr   = var(Pop_psth(1:ntr,:,:),[],3)./mean(Pop_psth(1:ntr,:,:),3);
%     
%     r_m = corrcoef(Pop_mean_tr');
%     r_v = corrcoef(Pop_var_tr');
%     r_f = corrcoef(Pop_ff_tr');
%     
%     figure; 
%     subplot(3,1,1);
%     histogram(r_m(:),-0.2:0.05:1.1)
%     hold on
%     plot(mean(r_m(:)),0,'og')
%     subplot(3,1,2)
%     histogram(r_v,-0.2:0.05:1.1)
%     hold on
%     plot(mean(r_v(:)),0,'og')
%     subplot(3,1,3)
%     histogram(r_f,-0.2:0.05:1.1)
%     hold on
%     plot(mean(r_f(:)),0,'og')
    
    
    %% Add stim, psth, and mahal dynamics to figure
    
    % Stimulus RMS
    subplot(hs(ist,1)); hold on
    plot(0:Duration, mean(stim,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','k')
    set(gca,'Color','none','xtick',[],'ytick',[])
    rmsrange = [rmsrange; get(gca,'ylim')];
    set(gca,'ylim',mode(rmsrange,1))
    
    % PSTH + Mahal
    subplot(hs(ist,2)); hold on
    plot(0:Duration, mean(psth,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','k')
    plot(0:Duration, mean(mdst,1,'omitnan'),...
        'LineWidth',psthlinewidth,'Color','m')
    set(gca,'Color','none','xtick',[],'ytick',[])
    ylim([0 ymaxval])
    
    
    %% Finish plot
    suptitle(sprintf('%s   |   %s %s   |   %i units',...
        Info.stim_ID_key{stid}, SUBJECT,SESSION,N_un ))
    
    subplot(hs(ist,3)); hold on
    set(gca,'ytick',add_y(2:end),'yticklabel',ytickstr,'ylim',[0 add_y(end)+1],...
        'xtick',unique([0:500:Duration Duration]),'xticklabel',unique([0:500:Duration Duration]))
    xlabel('Time (ms)')
    ylabel([num2str(mode(diff(add_y))) ' trials each'])    
%     ax = ancestor(gca,'axes');
%     yrule = ax.YAxis;
%     yrule.FontSize = 16;
    
    
    %% Save figure
    
    savedir = fullfile(fn.processed,'Sess_RastersMahal',sprintf('%s_%s',SUBJECT,SESSION));
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    savename = sprintf('%s_%s_%s',SUBJECT,SESSION,Info.stim_ID_key{ist});

    print_eps_kp(hf(ist),fullfile(savedir,savename))
    
%     set(gcf,'PaperOrientation','landscape')
%     print(hf(ist),fullfile(savedir,savename),'-dpdf','-bestfit')
    
    
end %ist





end