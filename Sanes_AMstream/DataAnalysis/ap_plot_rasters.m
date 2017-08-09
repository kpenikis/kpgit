function ap_plot_rasters(subject, session, channels, clus)
%
%  pp_plot_rasters(subject, session, channel, clu)
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%
%  KP, 2016-04; last updated 2017-06
%  (1,:) = Instantaneous AM rate <-- if Trials stim set, just this
%  (2,:) = Sound output          <-- if Trials stim set, just this
%  (3,:) = AM depth
%  (4,:) = dB SPL
%  (5,:) = HP
%  (6,:) = LP
%  (7,:) = Spout TTL
%  (8,:) = Block label
%


%!!!!!!!!!!!
SUonly = 0;
%!!!!!!!!!!!


% IF SAVING PDF FILES
figFontSize      = 14;
rasterMarkerSize = 10;
rasterLineWidth  = 0.5;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];
figsize2 = [1 scrsz(4) scrsz(3)/2 scrsz(4)];


histbinsize = 20;
anbinsize   = 10;


colors = hsv(6);
colors = [colors; 0.5.*hsv(4)];




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

% for irb = 7:10
%     find(cell2mat(cellfun(@eval,strsplit(Info.blockKey{irb}),'UniformOutput',0))==2)
% end


TotalDur_s = size(SoundData,2)/Info.fs_sound;


% Get unique dBSPLs
dBSPL = unique(SoundData(4,:));
rm_i=[];
for ii = 1:numel(dBSPL)
    if (numel(find(SoundData(4,:)==dBSPL(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
dBSPL(rm_i) = [];

% Get unique noisebands (based on LP)
LP = unique(SoundData(6,:));
rm_i=[];
for ii = 1:numel(LP)
    if (numel(find(SoundData(6,:)==LP(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
LP(rm_i) = [];

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
        close all
        
        try
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        % spiketrials = spikes.trials(unit_in);
        
        if isempty(spiketimes)
            error('no spike events found for this clu')
        elseif spikes.labels(spikes.labels(:,1)==clu,2) == 4
            warning('  this clu is labeled as noise. are you sure you want to plot?')
            keyboard
        end
        
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Set ymax based on overall FR of unit
        
        yvals = [50 100 200 300 400];
        yin = find( ( max(yvals, 2.5 * numel(spiketimes)/(size(SoundData,2)/Info.fs_sound)) - yvals )==0 );
        ymaxval = yvals(yin(1));
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        %%        
        
        % Also step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                for amd = AMdepth
                    
                    % Get this unit/stimulus combo's N clean blocks (trials) 
                   blocks_N = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    
                    %%
                    % Set up figure
                    fprintf('plotting ch %i clu %i...\n',channel,clu)
                    hf_b = figure;
                    set(hf_b,'Position',figsize1,'NextPlot','add')
                    hf_t = figure;
                    set(hf_t,'Position',figsize1,'NextPlot','add')
                    hf_ff = figure;
                    set(hf_ff,'Position',figsize1,'NextPlot','add')
                    
                    
                    % Create empty vectors for data
                    FR_pdc_pd = nan(90,4,6);
                    trtrR = nan(1,10);
                    FF = nan(1,10);
                    for ir = 1:numel(AMrates)
                        eval(sprintf( 'FRresp_IR_%i = [];',AMrates(ir) ))
                        eval(sprintf( 'FRresp_info_%i = [];',AMrates(ir) ))
                    end
                    
                    
                    % Go through each stiulus
                    for ib = 1:numel(Info.blockKey)
                        
                        %if unmodulated or silent, skip for now
                        if ib>=11, continue, end
                        
                        % Get samples of beginning and end of block
                        bkStart_samps = 1+find( diff(SoundData(8,:)==ib) ==  1 );
                        bkStop_samps  =   find( diff(SoundData(8,:)==ib) == -1 );
                        
                        % Remove blocks not this spl, lp noise, am depth,
                        % OR that contain artifact
                        rm_bk = [];
                        for it = 1:numel(bkStart_samps)
                            if ib==5 && (SoundData(8,bkStart_samps(it)-1)==11)
                                rm_bk = [rm_bk it];
                            elseif ib~=5 && (SoundData(8,bkStart_samps(it)-1)==11)
                                keyboard
                            end
                            if ~all(SoundData(4,bkStart_samps(it):bkStop_samps(it))==spl)...               
                                    || ~all(SoundData(6,bkStart_samps(it):bkStop_samps(it))==lpn)...      
                                    || ~all(SoundData(3,bkStart_samps(it):bkStop_samps(it))==amd)...      
                                    || any( intersect(bkStart_samps(it):bkStop_samps(it),ArtifactFlag))...  
                                    || ~all(SoundData(7,bkStart_samps(it):bkStop_samps(it))==1)          
                                rm_bk = [rm_bk it];
                            end
                        end
                        bkStart_samps(rm_bk) = [];
                        bkStop_samps(rm_bk) = [];
                        
                        % Convert to ms
                        bkStart_ms = round( bkStart_samps / Info.fs_sound*1000 );
                        bkStop_ms  = round( bkStop_samps  / Info.fs_sound*1000 );
                        
                        
                        if numel(bkStart_ms) ~= numel(bkStop_ms)
                            keyboard
                        end
                        %view preceeding blocks:  hpb = hist(SoundData(8,bkStart_samps-1),0:12)
                        
                        
%                         % Sort trials by preceding block
%                         [pbs,ipb]=sort(SoundData(8,bkStart_samps-1));
%                         bkStart_ms = bkStart_ms(ipb);
%                         bkStop_ms  = bkStop_ms(ipb);
%                         bkStart_samps = bkStart_samps(ipb);
%                         bkStop_samps  = bkStop_samps(ipb);
                        
                        

                        % Now prepare to collect response
                        
                        raster_x   = []; 
                        raster_y   = [];
                        stim       = nan(numel(bkStart_ms),max(bkStop_samps-bkStart_samps));
                        hist_raw   = zeros(numel(bkStart_ms),max(bkStop_ms-bkStart_ms));
                        FRresp     = nan(numel(bkStart_ms),1);
                        SpkTs_pdc  = [];
                        SpkTs_pdc2  = [];
                        
                        for it = 1:numel(bkStart_ms)
                            
                            % Sound rms envelope
                            stim(it,1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            
                            % Get overall spiking data for this block
                            sp=[]; sp = spiketimes( spiketimes>=bkStart_ms(it) & spiketimes<bkStop_ms(it) ) - bkStart_ms(it) +1;
                            hist_raw(it,sp) = 1;
                            raster_x = [raster_x sp];
                            raster_y = [raster_y it .* ones(1,numel(sp))];
                            % Save concatenated spiketimes for VS
%                             SpkTs_pdc2 = [SpkTs_pdc2 sp + (it-1)*(1000/AMrates(ib))];
                            
                            FRresp(it) = numel(sp) / ((bkStop_ms(it)-bkStart_ms(it))/1000);                            
                            
                            % Get spiking data for individual PERIODIC periods
                            if ib<=6 
                                % Estimate samples of period restarts
                                % (phase=0)
                                pd_starts = 0 : round(Info.fs_sound/AMrates(ib)) : (bkStop_samps(it)-bkStart_samps(it));
                                                                
                                % Get spikes
                                which_pd_sec = [0 0.5 1 1.5];
                                for ipd = 1:4
                                    [~,pdc_pd] = min(abs(pd_starts-which_pd_sec(ipd)*Info.fs_sound));
                                    pd_start_ms = round( ( pd_starts(pdc_pd) + bkStart_samps(it) ) /Info.fs_sound *1000 );
                                    FR_pdc_pd(it,ipd,ib) = sum( spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) / (1/AMrates(ib));
                                    SpkTs_pdc = [SpkTs_pdc spiketimes(spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) - pd_start_ms + (it-1)*(1000/AMrates(ib))];
                                end
                                
                                % Calculate VS, periodic blocks
                                [VS(ib),RS(ib),P(ib)] = vectorstrength( SpkTs_pdc, AMrates(ib) );
                                
                            end
                            
                            
                            % Get spiking data for each period of IR blocks
                            if ib>6 && ib<11
                                
                                find(cell2mat(cellfun(@eval,strsplit(Info.blockKey{ib}),'UniformOutput',0))==2);
                                
                                newRate = -1+bkStart_samps(it) + find(diff( [SoundData(1, (bkStart_samps(it)-1):(bkStop_samps(it)) ) 0] ));
                                if numel(newRate)~=13
                                    keyboard
                                end
                                for ir = 1:(numel(newRate)-1)
                                    this_rate = SoundData(1, newRate(ir));
                                    t = round( [newRate(ir) newRate(ir+1)-1] / Info.fs_sound*1000 );
                                    sp=[]; sp = spiketimes( spiketimes>=t(1) & spiketimes<=t(2) ) - t(1) +1;
                                    eval(sprintf('FRresp_IR_%i = [FRresp_IR_%i; numel(sp) / (diff(t)/1000)];',this_rate, this_rate ))
                                    eval(sprintf('FRresp_info_%i = [FRresp_info_%i; (newRate(ir)-bkStart_samps(it))/Info.fs_sound ib];', this_rate, this_rate ))
                                end
                                
                                % Calculate VS, IR blocks
                                
                            end
                            
                            
                            
                        end %it
                        
                        
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        % Calculate trial-trial correlations
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        % First, bin trial spike counts
                        hist_bin = binspikecounts(hist_raw,anbinsize);
                        
                        % Now choose N random pairs of trials from hist_bin
                        % and average their correlation coefficients
                            % N is determined by stim with fewest trials
                        bk_i = randperm(size(hist_bin,1));
                        fullR = corrcoef(hist_bin(bk_i(1:min(blocks_N)),:));
                        trtrR(ib) = mean(fullR(fullR~=1),'omitnan');
                        
                        
                        % And the average FF across bins
                        FF(ib) = mean(var(hist_bin,1)./mean(hist_bin,1),'omitnan');
                        FFsem(ib) = std(var(hist_bin,1)./mean(hist_bin,1),'omitnan')/sqrt(size(hist_bin,2));
                        
                        
                        
                        %% Plots
                        
                        if ib<7
                            blockstr = [Info.blockKey{ib} 'Hz'];
                        else
                            blockstr = Info.blockKey{ib};
                        end
                        unType = {'' 'SU' 'MU'};
                        titlestr1 = sprintf('Block %s\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz | depth=%i',...
                            blockstr,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn,amd);
                        titlestr2 = sprintf('%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz | depth=%i',...
                            subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn,amd);
                        titlestr3 = sprintf('FR individual periods, by context and time in block\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz | depth=%i',...
                            subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn,amd);
                        
                        
                        % RASTER/PSTH
                        
                        figure;
                        set(gcf,'Position',figsize1)
                        subplot(5,1,1:2); hold on
                        plot(  raster_x  ,  raster_y  , 'k+','MarkerSize',rasterMarkerSize, 'LineWidth', rasterLineWidth)
                        for ii=1:it
                            plot([0 max(bkStop_ms-bkStart_ms)], [ii ii], 'k', 'LineWidth', rasterLineWidth)
                        end
                        ylim([0 it+1])
                        xlim([0 max(bkStop_ms-bkStart_ms)])
                        set(gca,'xtick',[])
                        ylabel('trials')
                        set(gca,'FontSize',figFontSize)
                        
                        title(titlestr1)
                        
                        subplot(5,1,3)
                        plot(stim(:,1:round(Info.fs_sound/1000):end)','Color',colors(ib,:))
                        xlim([0 size(stim(:,1:round(Info.fs_sound/1000):end),2)])
                        set(gca,'xtick',[],'ytick',[])
                        ylabel('stim rms')
                        set(gca,'FontSize',figFontSize)
                        
                        
                        subplot(5,1,4:5); hold on
                        hist_smooth = binspikecounts(hist_raw,histbinsize)./histbinsize*1000;
                        plot(histbinsize:histbinsize:size(hist_raw,2),mean(hist_smooth,1),'k','LineWidth',3)
                        ylim([0 ymaxval])
                        xlabel('time (ms)')
                        ylabel('spikes/s')
                        set(gca,'FontSize',figFontSize)
                        
                        
                        % SAVE RASTER/PSTH FIGURE
                        
                        savedir = fullfile(fn.rasterplots,['ch' num2str(channel)]);
                        if ~exist(savedir,'dir')
                            mkdir(savedir)
                        end
                        savename = sprintf('%s_%s_ch%i_clu%i_%s_bk%s_%idB_LP%ihz',...
                            subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},blockstr,spl,lpn);
                        set(gcf,'PaperOrientation','landscape');
                        print(gcf,'-dpdf',fullfile(savedir,savename),'-bestfit')
                        
                        
                        
                        %% Plot FR response for this block type
                        
                        figure(hf_b); 
                        hold on
                        ee=errorbar(ib,mean(FRresp,'omitnan'),std(FRresp,'omitnan')/sqrt(numel(FRresp)),...
                            'o','MarkerSize',16,'LineWidth',2,...
                            'Color',colors(ib,:),'MarkerFaceColor',colors(ib,:),'MarkerEdgeColor','none');
                        
                        if ib>6
                            ee.MarkerEdgeColor = [0 0 0];
                        end
                        
                        meanFRs(ib) = mean(FRresp,'omitnan');
                        
                        
                        % Plot average trial-trial correlation
                        figure(hf_t);
                        hold on
                        ee=plot(ib,trtrR(ib),...
                            'o','MarkerSize',16,'LineWidth',2,...
                            'Color',colors(ib,:),'MarkerFaceColor',colors(ib,:),'MarkerEdgeColor','none');
                        
                        if ib>6
                            ee.MarkerEdgeColor = [0 0 0];
                        end
                        
                        % Plot FF and spread across bins
                        figure(hf_ff);
                        hold on
                        ee=errorbar(ib,FF(ib),FFsem(ib),...
                            'o','MarkerSize',16,'LineWidth',2,...
                            'Color',colors(ib,:),'MarkerFaceColor',colors(ib,:),'MarkerEdgeColor','none');
                        
                        if ib>6
                            ee.MarkerEdgeColor = [0 0 0];
                        end
                        
                        
                        
                    end %ib
                    
                    % Finish FR fig
                    figure(hf_b)
                    % Plot mean of means for PDC and IR
                    plot([0.5 6.5], [mean(meanFRs(1:6))  mean(meanFRs(1:6))] ,'k--')
                    plot([6.5 10.5],[mean(meanFRs(7:10)) mean(meanFRs(7:10))],'k--')
                    
                    % Add some formatting 
                    title(sprintf('FR MTF\n%s',titlestr2))
                    ylim([0 0.6*ymaxval])
                    xlim([0 13])
                    ylabel('mean FR')
                    set(gca,'xtick',1:numel(Info.blockKey))
                    xticklabels(Info.blockKey)
                    xtickangle(45)
                    set(gca,'TickLabelInterpreter','none')
                    set(gca,'FontSize',figFontSize)
                    
                    % Finish tr-tr corr fig
                    figure(hf_t)
                    % Plot mean of means for PDC and IR
                    plot([0.5 6.5], [mean(trtrR(1:6))  mean(trtrR(1:6))] ,'k--')
                    plot([6.5 10.5],[mean(trtrR(7:10)) mean(trtrR(7:10))],'k--')
                    
                    % Add some formatting
                    title(sprintf('Trial Corr MTF\n%s',titlestr2))
                    ylims = get(gca,'ylim');
                    ylim([0 ylims(2)*1.2])
                    xlim([0 13])
                    ylabel(['Average correlation coeff of single trial pairs, ' num2str(anbinsize) ' ms bins'])
                    set(gca,'xtick',1:numel(Info.blockKey))
                    xticklabels(Info.blockKey)
                    xtickangle(45)
                    set(gca,'TickLabelInterpreter','none')
                    set(gca,'FontSize',figFontSize)
                    
                    % Finish tr-tr corr fig
                    figure(hf_ff)
                    plot([0 13], [1 1] ,'Color',[0.5 0.5 0.5])
                    % Plot mean of means for PDC and IR
                    plot([0.5 6.5], [mean(FF(1:6))  mean(FF(1:6))] ,'k--')
                    plot([6.5 10.5],[mean(FF(7:10)) mean(FF(7:10))],'k--')
                    
                    % Add some formatting
                    title(sprintf('FF MTF\n%s',titlestr2))
                    ylims = get(gca,'ylim');
                    ylim([1-max(abs(ylims-1)) 1+max(abs(ylims-1))])
                    xlim([0 13])
                    ylabel(['FF, average across ' num2str(anbinsize) ' ms bins'])
                    set(gca,'xtick',1:numel(Info.blockKey))
                    xticklabels(Info.blockKey)
                    xtickangle(45)
                    set(gca,'TickLabelInterpreter','none')
                    set(gca,'FontSize',figFontSize)
                    
                    

                    % SAVE MTF FIGURES
                    savedir = fn.anplots;
                    if ~exist(savedir,'dir')
                        mkdir(savedir)
                    end
                    % FR
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_MTF-FR_%idB_LP%ihz',...
                        subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn);
                    set(hf_b,'PaperOrientation','landscape');
                    print(hf_b,'-dpdf',fullfile(savedir,savename),'-bestfit')
                    % tr-tr corr
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_MTF-trCorr_%idB_LP%ihz',...
                        subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn);
                    set(hf_t,'PaperOrientation','landscape');
                    print(hf_t,'-dpdf',fullfile(savedir,savename),'-bestfit')
                    % FF
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_MTF-FF_%idB_LP%ihz',...
                        subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn);
                    set(hf_ff,'PaperOrientation','landscape');
                    print(hf_ff,'-dpdf',fullfile(savedir,savename),'-bestfit')
                    
                    
                    
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    % Barplots of individual AM rates from the IR blocks
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    hf_p = figure;
                    clf
                    set(gcf,'NextPlot','add','Position',figsize2)
                    
                    for ir = AMrates
                        
                        subplot(3,2,find(ir==AMrates)); hold on
                        
                        % Periodic blocks
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
                        errorbar([0 0.5 1 1.5],mean(FR_pdc_pd(:,:,ir==AMrates),1,'omitnan'),std(FR_pdc_pd(:,:,ir==AMrates),1,'omitnan')/sqrt(sum(~isnan(FR_pdc_pd(:,1,ir==AMrates)))),...
                            'd-','MarkerSize',16,'LineWidth',2,...
                            'Color',colors(ir==AMrates,:),'MarkerFaceColor',colors(ir==AMrates,:),'MarkerEdgeColor','none');
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
                        
                        % IR blocks
                        IRdata = eval(sprintf('FRresp_IR_%i',ir));
                        IRinfo = eval(sprintf('FRresp_info_%i',ir));
                        
                        % Separate by which IR block
                        for ib=7:10
                            
                            % and by early/late pds
                            pdStarts = unique(sign( IRinfo(IRinfo(:,2)==ib,1) - sum(1./AMrates) ))';
                            if numel(pdStarts)~=2, keyboard, end
                            
                            for ipd = pdStarts
                                
                                plotdata = IRdata( sign(IRinfo(:,1)-sum(1./AMrates))==ipd & IRinfo(:,2)==ib );
                                if isempty(plotdata), keyboard, end
                                
                                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
                                errorbar(mean(IRinfo( sign(IRinfo(:,1)-sum(1./AMrates))==ipd & IRinfo(:,2)==ib ,1)),...
                                    mean(plotdata), std(plotdata)/sqrt(numel(plotdata)),...
                                    'd','MarkerSize',16,'LineWidth',2,...
                                    'Color',colors(ib,:),'MarkerFaceColor',colors(ib,:),'MarkerEdgeColor','none');
                                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
                                
                            end %ipd
                        end %ib
                        
                        xlim([0 2])
                        ylim([0 ymaxval])
                        if ir==32
                            ylabel('mean FR for individual periods')
                            xlabel('pd start time within block (s)')
                            set(gca,'FontSize',figFontSize+2)
                        else
                            set(gca,'xticklabel',[],'yticklabel',[])
                        end
                        
                    end %ir
                    
                    % Add super title
                    suptitle(titlestr3)
                    
                    
                    % SAVE FIGURE
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_pdFR-time_%idB_LP%ihz',...
                        subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn);
                    %     print(hf_b,fullfile(savedir,savename),'-depsc','-tiff')
                    set(hf_p,'PaperOrientation','landscape');
                    print(hf_p,'-dpdf',fullfile(savedir,savename),'-bestfit')
                    
                    
                    
                    pause(2)
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end %function




