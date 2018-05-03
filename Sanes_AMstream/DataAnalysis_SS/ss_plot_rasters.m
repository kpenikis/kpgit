function ss_plot_rasters(subject, session, channels, clus)
%
%  pp_plot_rasters(subject, session, channel, clu)
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%
%  KP, 2017-08
%  (1,:) = Instantaneous AM rate <-- if Trials stim set, just this
%  (2,:) = Sound output          <-- if Trials stim set, just this
%  (3,:) = AM depth
%  (4,:) = dB SPL
%  (5,:) = CF
%  (6,:) = Spout TTL
%  (7,:) =   x
%  (8,:) = Block label
%


%!!!!!!!!!!!!!!!!!
SUonly   =  0;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10; 
%!!!!!!!!!!!!!!!!!


% IF SAVING PDF FILES
figFontSize      = 14;
rasterMarkerSize = 50;
rasterLineWidth  = 0.5;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];

hf_bar = figure;
set(hf_bar,'Position',figsize)


histbinsize = 20;
anbinsize   = 10;
smthbinsize = 50;


col_0  = [ 16  20 100]./255;
col_1  = [120  16  16]./255;
col_01 = [225 104  26]./255;
col_10 = [ 30 122 245]./255;



%%
% Load data files
fn = set_paths_directories(subject,session);
fprintf('loading data...\n')
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_SoundData',subject,session); load(fullfile(fn.processed,subject,filename));
load(fullfile(fn.stim,'IRsequences.mat'))


%%
% GET STIM INFO

AMrates = [2 4 8 16 32 64];
fs = round(Info.fs_sound);

% Get unique dBSPLs
dBSPL = unique(SoundData(4,:));
rm_i=[];
for ii = 1:numel(dBSPL)
    if (numel(find(SoundData(4,:)==dBSPL(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
dBSPL(rm_i) = [];

% Get unique AM depths
AMdepth = unique(SoundData(3,:));
rm_i=[];
for ii = 1:numel(AMdepth)
    if (numel(find(SoundData(3,:)==AMdepth(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
AMdepth(rm_i) = [];





%% GET SPIKE TIMES

% Step through all channels if not specified
if nargin<3 && ~exist('channels','var')
    channels = [1:7 9:16];
end


% STEP THROUGH EACH CHANNEL
for channel = channels
    
    % Artifact for this channel
    ArtifactFlag = Info.artifact(channel).SDsamples;
    
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if all(spikes.labels(:,2)==1)
            disp(' SESSION MAY NOT BE MANUALLY SORTED YET')
            return
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
            disp(' !! no valid clus for this channel')
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
        end
    end
    
    % STEP THROUGH EACH CLU
    
    for clu = clus'
        
        % !! Only SU for now !!
        if SUonly && (spikes.labels(spikes.labels(:,1)==clu,2) ~= 2)
            continue
        end
        
        try
            spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        
        if isempty(spiketimes)
            error('no spike events found for this clu')
        elseif spikes.labels(spikes.labels(:,1)==clu,2) == 4
            warning('  this clu is labeled as noise. are you sure you want to plot?')
            keyboard
        end
        
        fprintf('plotting ch %i clu %i...\n',channel,clu)
        if numel(spiketimes) < round(FRcutoff*length(SoundData)/Info.fs_sound)
            disp(' few spiking events')
        end
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Make constant stream of 
        Stream_Spks = zeros(1,1000*ceil((size(SoundData,2)/Info.fs_sound)));
        Stream_Spks(spiketimes) = 1;
        
        %either with standard 20 ms bin
%         Stream_FRbin = 1000*(binspikecounts(Stream_Spks,histbinsize)/histbinsize);
%         Stream_FRbin(isinf(Stream_FRbin)) = nan;
%         foo = repmat(Stream_FRbin,histbinsize,1);
%         Stream_FR = reshape(foo,1,histbinsize*length(Stream_FRbin));
%         Stream_FR = Stream_FR(1:ceil(length(SoundData)/Info.fs_sound*1000));
        
        %or with sliding 50 ms boxcar
        Stream_FRsmooth = smoothFR(Stream_Spks,smthbinsize);
        Stream_FRsmooth = Stream_FRsmooth(1:ceil(length(SoundData)/Info.fs_sound*1000));
        
        % Determine the time of the beginning of actual recording session
        % max of first spiketime or 5 seconds before unmodulated sound came
        % on
        sampStart = find(diff(SoundData(8,:))==-1);
        msStart   = max( spiketimes(1), round(sampStart(1)/Info.fs_sound*1000)-5000 );
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Convert FR to z-score
        
        REFERENCE = 'unmod';'stream';
        markRED = 0;
        
        switch REFERENCE
            
            case 'stream'
                Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
            
            case 'unmod'
                
%                 if ~any(strcmp(Zdata.Properties.VariableNames,'preUNM'))
%                     Zdata.preUNM   = 0;
%                     Zdata.postUNM  = 0;
%                 end
                
                % Get samples of unmodulated sound and of silence
                unmod_samps  = find(SoundData(8,:)==11);
                
                % Split unmod portion into before and after
                try
                    if ~isempty(find(diff(unmod_samps)>1))
                        unmod1_ms = round( [unmod_samps(1)   unmod_samps(0+find(diff(unmod_samps)>1))] /Info.fs_sound*1000 );
                        unmod2_ms = round( [unmod_samps(1+find(diff(unmod_samps)>1)) unmod_samps(end)] /Info.fs_sound*1000 );
                        
%                         diff(unmod1_ms)/1000
%                         diff(unmod2_ms)/1000
                        
                        meanFR = mean([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        stdFR = std([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        
                    else
                        unmod1_ms = round( [unmod_samps(1) unmod_samps(end)]/Info.fs_sound*1000 );
                        unmod2_ms = [0 0];
                        
                        meanFR = mean(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        stdFR = std(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        
                        fprintf(' !! not much unmod data for sess %s ch%i clu%i \n',session,channel,clu)
%                         markRED = 1;
                    end
                    
                catch
                    keyboard
                end
                
                Stream_zscore = (Stream_FRsmooth(msStart:end) - meanFR) / stdFR;
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
                
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Set ymax based on overall FR of unit
        
        yvals = [50 100 200 300 400];
        yin = find( ( max(yvals, 2.5 * numel(spiketimes)/(size(SoundData,2)/Info.fs_sound)) - yvals )==0 );
        ymaxval = yvals(yin(1));
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        %%  SOUND PARAMS
        
        % Also step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for amd = AMdepth
                
                % Get this unit/stimulus combo's N clean blocks (trials)
                [blocks_N,minTime] = get_N_clean_blocks_SS(SoundData,Info,ArtifactFlag,spl,amd);
                
                
                
                %%  BLOCKS
                
                legx = nan(numel(find(blocks_N)),1);
                legstr = cell(numel(find(blocks_N)),1);
                nb = 1;
                
                % Go through each stimulus
                for ib = find(blocks_N)
                    
                    
                    % Skip if few trials
                    if blocks_N(ib)<minTrs
                        continue
                    end
                    
                    % Get this block start times
                    [bkStart_samps,bkStop_samps,...
                        bkStart_ms,bkStop_ms]  =  get_blockOnsets_SS( SoundData,...
                        ib,spl,amd,ArtifactFlag,Info.fs_sound);
                    
                    
                    % Now prepare to collect response
                    
                    raster_x_0     = [];
                    raster_y_0     = [];
                    stim_0         = nan(1,max(bkStop_samps-bkStart_samps));
                    psth_smooth_0  = nan(1,max(bkStop_ms-bkStart_ms));
                    trs_0          = [];
                    
                    raster_x_1     = [];
                    raster_y_1     = [];
                    stim_1         = nan(1,max(bkStop_samps-bkStart_samps));
                    psth_smooth_1  = nan(1,max(bkStop_ms-bkStart_ms));
                    trs_1          = [];
                    
                    raster_x_01    = [];
                    raster_y_01    = [];
                    stim_01        = nan(1,max(bkStop_samps-bkStart_samps));
                    psth_smooth_01 = nan(1,max(bkStop_ms-bkStart_ms));
                    sampTrans_01   = [];
                    trs_01         = [];
                    
                    raster_x_10    = [];
                    raster_y_10    = [];
                    stim_10        = nan(1,max(bkStop_samps-bkStart_samps));
                    psth_smooth_10 = nan(1,max(bkStop_ms-bkStart_ms));
                    sampTrans_10   = [];
                    trs_10         = [];
                    
                    stim           = nan(numel(bkStart_ms),ceil(minTime/1000*fs));
                    
                    
                    % Get data from each trial of this block type
                    
                    for it = 1:numel(bkStart_ms)
                        
                        stim(it,:) = envelope(SoundData(2, bkStart_samps(it)+[1:ceil(minTime/1000*fs)]),20,'rms');
                        
                        % Get spiketimes for this block
                        sp=[]; sp = spiketimes( spiketimes>=bkStart_ms(it) & spiketimes<bkStop_ms(it) ) - bkStart_ms(it) +1;
                        
                        % Decide where to collect this data, based on what
                        % happened with the CF 
                        
                        if sum(diff( SoundData(5,bkStart_samps(it):bkStop_samps(it)) ))==1
                            % CF  0 to 1
                            trs_01(end+1) = it;
                            sampTrans_01(end+1) = find(diff( SoundData(5,bkStart_samps(it):bkStop_samps(it)) ) == 1);
                            
                            stim_01(numel(trs_01),1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            raster_x_01 = [raster_x_01 sp];
                            raster_y_01 = [raster_y_01 numel(trs_01) .* ones(1,numel(sp))];
                            psth_smooth_01(numel(trs_01),1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_FRsmooth(bkStart_ms(it):bkStop_ms(it)-1);
                            
                        elseif sum(diff( SoundData(5,bkStart_samps(it):bkStop_samps(it)) ))==-1
                            % CF  1 to 0
                            trs_10(end+1) = it;
                            sampTrans_10(end+1) = find(diff( SoundData(5,bkStart_samps(it):bkStop_samps(it)) ) == -1);
                            
                            stim_10(numel(trs_10),1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            raster_x_10 = [raster_x_10 sp];
                            raster_y_10 = [raster_y_10 numel(trs_10) .* ones(1,numel(sp))];
                            psth_smooth_10(numel(trs_10),1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_FRsmooth(bkStart_ms(it):bkStop_ms(it)-1);
                            
                        elseif sum(SoundData(5,bkStart_samps(it):bkStop_samps(it))) == 0
                            % CF just 0
                            trs_0(end+1) = it;
                            
                            stim_0(numel(trs_0),1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            raster_x_0 = [raster_x_0 sp];
                            raster_y_0 = [raster_y_0 numel(trs_0) .* ones(1,numel(sp))];
                            psth_smooth_0(numel(trs_0),1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_FRsmooth(bkStart_ms(it):bkStop_ms(it)-1);
                            
                        elseif sum(SoundData(5,bkStart_samps(it):bkStop_samps(it))) > 500 && sum(diff( SoundData(5,bkStart_samps(it):bkStop_samps(it)) ))==0
                            % CF just 1
                            trs_1(end+1) = it;
                            
                            stim_1(numel(trs_1),1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            raster_x_1 = [raster_x_1 sp];
                            raster_y_1 = [raster_y_1 numel(trs_1) .* ones(1,numel(sp))];
                            psth_smooth_1(numel(trs_1),1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_FRsmooth(bkStart_ms(it):bkStop_ms(it)-1);
                            
                        end
                        
                        
                    end %it
                    
                    if range(sampTrans_01)>2 || range(sampTrans_10)>2
                        keyboard
                    end
                    
                    
                    %% PLOT
                    
                    if ib<7
                        blockstr = [Info.blockKey{ib} 'Hz'];
                    else
                        blockstr = Info.blockKey{ib};
                    end
                    unType = {'' 'SU' 'MU'};
                    
                    titlestr = sprintf('Block %s (%i dBSPL)\n%s | %s | ch %i, clu %i (%s)',...
                        blockstr,spl,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)});
                    
                    
                    % Plot figure
                    
                    figure;
                    set(gcf,'Position',figsize)
                    
                    % RASTER
                    subplot(5,1,1:2); hold on
                    trcount = 0;
                    scatter(  raster_x_0  ,  raster_y_0  , rasterMarkerSize,...
                        'o', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor', col_0, 'MarkerEdgeColor','none')
                    trcount = trcount + max([raster_y_0 0]);
                    
                    scatter(  raster_x_1  ,  raster_y_1   + trcount  , rasterMarkerSize,...
                        'o', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor', col_1, 'MarkerEdgeColor','none')
                    trcount = trcount + max([raster_y_1 0]);
                    
                    scatter(  raster_x_01  ,  raster_y_01 + trcount , rasterMarkerSize,...
                        'o', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor', col_01, 'MarkerEdgeColor','none')
                    trcount = trcount + max([raster_y_01 0]);
                    
                    scatter(  raster_x_10  ,  raster_y_10 + trcount , rasterMarkerSize,...
                        'o', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor', col_10, 'MarkerEdgeColor','none')
                    trcount = trcount + max([raster_y_10 0]);
%                     for ii=1:it
%                         plot([0 max(bkStop_ms-bkStart_ms)], [ii ii], 'k', 'LineWidth', rasterLineWidth)
%                     end
                    plot([mean([sampTrans_10 sampTrans_01]) mean([sampTrans_10 sampTrans_01])]/Info.fs_sound*1000, [0 trcount+1], 'k', 'LineWidth', rasterLineWidth)
                    ylim([0 trcount+1])
                    xlim([0 minTime])
                    set(gca,'xtick',[])
                    ylabel('trials')
                    set(gca,'FontSize',figFontSize)
                    
                    title(titlestr)
                    
                    subplot(5,1,3); hold on
%                     plot(stim(:,1:round(fs/1000):end)','Color','k')
                    plot(stim_0','Color',col_0)
                    plot(stim_1','Color',col_1)
                    plot(stim_01','Color',col_01)
                    plot(stim_10','Color',col_10)
                    
                    xlim([0 size(stim_0,2)])
                    set(gca,'xtick',[],'ytick',[])
                    ylabel('stim rms')
                    set(gca,'FontSize',figFontSize)
                    
                    
                    % PSTHs
                    subplot(5,1,4:5); hold on
                    plot([mean([sampTrans_10 sampTrans_01]) mean([sampTrans_10 sampTrans_01])]/fs*1000, [0 ymaxval], 'k', 'LineWidth', rasterLineWidth)
                    ip(1)=plot(mean(psth_smooth_0, 1,'omitnan'),'Color',col_0,'LineWidth',3);
                    ip(2)=plot(mean(psth_smooth_1, 1,'omitnan'),'Color',col_1,'LineWidth',3);
                    ip(3)=plot(mean(psth_smooth_01,1,'omitnan'),'Color',col_01,'LineWidth',3);
                    ip(4)=plot(mean(psth_smooth_10,1,'omitnan'),'Color',col_10,'LineWidth',3);
                    ylim([0 ymaxval])
                    xlim([0 minTime])
                    xlabel('time (ms)')
                    ylabel('spikes/s')
                    set(gca,'FontSize',figFontSize)
                    legend(ip,{'low' 'high' 'low-high' 'high-low'},'Location','best')
                    
                    
                    % SAVE RASTER/PSTH FIGURE
                    
                    savedir = fullfile(fn.rasterplots,['ch' num2str(channel)]);
                    if ~exist(savedir,'dir')
                        mkdir(savedir)
                    end
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_bk%s_%idB',...
                        subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},blockstr,spl);
                    print_eps_kp(gcf,fullfile(savedir,savename),1);
%                     set(gcf,'PaperOrientation','landscape');
%                     print(gcf,'-dpdf',fullfile(savedir,savename),'-bestfit')
                    
                    
                    
                    %% QUANTIFY
                    
                    % Calculate avg FR response during first 250 ms after
                    % CF switch time
                    ms_win = 249;
                    transTime = round(mean([sampTrans_10 sampTrans_01])/fs*1000);
                    
                    FR_0     = mean( mean(psth_smooth_0(:, transTime+[0:ms_win]), 1,'omitnan'));
                    FRsem_0  = std( mean(psth_smooth_0(:,  transTime+[0:ms_win]), 2,'omitnan'))/sqrt(size(psth_smooth_0,1));
                    FR_1     = mean( mean(psth_smooth_1(:, transTime+[0:ms_win]), 1,'omitnan'));
                    FRsem_1  = std( mean(psth_smooth_1(:,  transTime+[0:ms_win]), 2,'omitnan'))/sqrt(size(psth_smooth_1,1));
                    FR_01    = mean( mean(psth_smooth_01(:,transTime+[0:ms_win]), 1,'omitnan'));
                    FRsem_01 = std( mean(psth_smooth_01(:, transTime+[0:ms_win]), 2,'omitnan'))/sqrt(size(psth_smooth_01,1));
                    FR_10    = mean( mean(psth_smooth_10(:,transTime+[0:ms_win]), 1,'omitnan'));
                    FRsem_10 = std( mean(psth_smooth_10(:, transTime+[0:ms_win]), 2,'omitnan'))/sqrt(size(psth_smooth_10,1));
                    
                    
                    figure(hf_bar); hold on
                    plot([nb+0 nb+1], [FR_0 FR_10],'Color',col_0,'LineWidth',2)
                    plot([nb+2 nb+3], [FR_1 FR_01],'Color',col_1,'LineWidth',2)
                    bp(1)=errorbar(nb+0, FR_0, FRsem_0,'s','Color',col_0,'LineWidth',2,...
                        'MarkerSize',30,'MarkerFaceColor',col_0,'MarkerEdgeColor',col_0);
                    bp(2)=errorbar(nb+1, FR_10, FRsem_10,'s','Color',col_10,'LineWidth',2,...
                        'MarkerSize',30,'MarkerFaceColor',col_10,'MarkerEdgeColor',col_10);
                    bp(3)=errorbar(nb+2, FR_1, FRsem_1,'s','Color',col_1,'LineWidth',2,...
                        'MarkerSize',30,'MarkerFaceColor',col_1,'MarkerEdgeColor',col_1);
                    bp(4)=errorbar(nb+3, FR_01, FRsem_01,'s','Color',col_01,'LineWidth',2,...
                        'MarkerSize',30,'MarkerFaceColor',col_01,'MarkerEdgeColor',col_01);
                    
                    legx(ib==find(blocks_N))   = nb+1.5;
                    legstr{ib==find(blocks_N)} = blockstr;
                    
                    nb = nb+6; 
                    
                    
                end %ib
                
                set(gca,'xtick',legx,'xticklabel',legstr)
                ylabel('average FR (Hz)')
                titlestr = sprintf('avg FR in %i ms after CF shift\n%s | %s | ch %i, clu %i (%s)',...
                    ms_win+1,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)});
                title(titlestr)
                legend(bp,{'low' 'high-low' 'high' 'low-high'},'Location','northeast','Orientation','horizontal')
                
                
                % Save barplot
                savedir = fullfile(fn.processed,'SS');
                if ~exist(savedir,'dir')
                    mkdir(savedir)
                end
                savename = sprintf('%s_%s_ch%i_clu%i_%s_FRbar%ims',...
                    subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},ms_win+1);
                print_eps_kp(gcf,fullfile(savedir,savename),1);
                
                
                
            end %amd
        end %spl
        
    end %clu
end %channel

end %function




