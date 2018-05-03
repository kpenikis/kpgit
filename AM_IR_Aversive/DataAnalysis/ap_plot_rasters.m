function ap_plot_rasters(subject, session, channels, clus)
%
%  pp_plot_rasters(subject, session, [channel, clu] )
%    Plots a raster and psth for each stimulus, separating trials by
%    the preceding stimulus. 
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2016-04; last updated 2018-03
%


%!!!!!!!!!!!
SUonly = 1;
%!!!!!!!!!!!
CalcZ  = 0;
nIterations = 100;
PLOT_RND = 0;
%!!!!!!!!!!!
PLOT_ALL = 1;
%!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
PLOT_FF = 0;
%!!!!!!!!!!!!!!!!!

rng('shuffle')
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
figsize1 = [1 2*scrsz(4)/3 scrsz(3) 2*scrsz(4)/3];
figsize2 = [1 1 scrsz(3) scrsz(4)/3];

colors = [  0   0   0;...
           84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [colors; 0.7.*bone(2)];
% previous expt : 64 Hz (magenta) : 255 205  60 


%%
% Load data files
fn = set_paths_directories(subject,session,1);
fprintf('loading data...\n')
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(fn.processed,subject,filename));

% KS
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],24414)
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],Info.fs)


% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);


%% STEP THROUGH CHANNELS AND CLUS

% Step through all channels if not specified
if nargin<3 && ~exist('channels','var')
    channels = [1:7 9:16];
end

for channel = channels
        
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
%             disp(' !! no valid clus for this channel')
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
        end
    end
    
    % STEP THROUGH EACH CLU
    
    for clu = clus'
        
        close all
        
        % GET SPIKETIMES
        
        if SUonly && (spikes.labels(spikes.labels(:,1)==clu,2) ~= 2)
            continue
        end
        
        switch spikes.labels(spikes.labels(:,1)==clu,2)
            case 2
                unType = 'SU';
            case 3
                unType = 'MU';
        end
        
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        % spiketrials = spikes.trials(unit_in);
        
        if isempty(spiketimes)
            error('no spike events found for this clu')
        elseif spikes.labels(spikes.labels(:,1)==clu,2) == 4
            warning('  this clu is labeled as noise. are you sure you want to plot?')
            keyboard
        end
        fprintf('plotting ch %i clu %i...\n',channel,clu)
%         
%         figure;
%         hist(diff(spiketimes(spiketimes>TrialData.onset(1))),0:1:257)
%         xlim([0 256])
        
        
        %% Plot of smoothed FR over entire session
        
%         bs_smth = 5000;
%         Stream_FRsmooth = convertSpiketimesToFR(spiketimes,...
%             length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
%         
%         figsize2 = [1 scrsz(4)/3 scrsz(3) scrsz(4)/3];
%         figure;
%         set(gcf,'Position',figsize2,'NextPlot','add')
%         plot((1:length(Stream_FRsmooth))./1000,Stream_FRsmooth,'k')
%         hold on
%         ylimits=get(gca,'ylim');
%         plot([TrialData.onset(2) TrialData.onset(2)]./1000,ylimits,':b','LineWidth',2)
%         fill([TrialData.onset(1) TrialData.onset(1) TrialData.offset(1) TrialData.offset(1)]./1000,[0 1 1 0].*ylimits(2),...
%             'b','FaceAlpha',0.3,'EdgeColor','none')
%         plot((1:length(Stream_FRsmooth))./1000,Stream_FRsmooth,'k')
%         xlim([0 length(Stream_FRsmooth)./1000])
%         xlabel('Time in session (seconds)')
%         ylabel('FR (spikes/s)')
        
%         fill([1:length(SpoutStream) length(SpoutStream):-1:1]./1000,[SpoutStream zeros(size(SpoutStream))].*ylimits(2),...
%             'g','FaceAlpha',0.3,'EdgeColor','none')
%         
        

        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Check the duration of silence at the beginning
        if ((TrialData.offset(1) - TrialData.onset(1))/1000) < 15
            keyboard
        end
        
        bs_smth = 20;
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        
        
        % Check if unit has very low overall firing rate
        ExcludeUnit = 0;
        if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
            disp(' few spiking events')
            ExcludeUnit = 1;
        end
        
        
        
        %%
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
                [all_TDidx,Ntrials1,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                
                allStim = unique(TrialData.trID(all_TDidx));
                meanFR_all = nan(1,numel(allStim));
                stimTr_nSpk_minD = nan(max(Ntrials1),numel(allStim));
                stimTr_FR_Norm   = nan(max(Ntrials1),numel(allStim));
                eachTr_nSpk_minD = []; 
                eachTr_nSpk_Norm = []; 
                eachTr_stid_iti  = []; 
                eachTr_pstid     = []; 
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                for stid = 1:numel(allStim)
                    
                    
                    %% FIRST, COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==allStim(stid));
                    
                    %%%  plot for Trial and ITI stimuli separately
                    ITIflag = unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for is = 1:numel(ITIflag)
                        
                        st_TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
                        
                        pst_TDidx = nan(size(st_TDidx));
                        skip_it = [];
                        
                        for it = 1:numel(st_TDidx)
                            
                            % Get TD index of previous (clean) trial
                            
                            if (find(all_TDidx==st_TDidx(it))-1)==0
                                skip_it = [skip_it it];
                                continue
                            end
                            pst_TDidx(it) = all_TDidx(find(all_TDidx==st_TDidx(it))-1);
                            
                            % Check how long between the onset of this trial
                            % and the offset of the previous trial. If there is
                            % a gap between them, skip this transition.
                            if (TrialData(st_TDidx(it),:).onset - TrialData(pst_TDidx(it),:).offset) > 1%ms
                                skip_it = [skip_it it];
                            end
                        end %it
                        
                        % Now have final trials to work with
                        st_TDidx(skip_it)  = [];
                        pst_TDidx = st_TDidx-1;
                        
                        % Get timestamps and durations
                        clear t1 t2 t3 Durations t_win 
                        t2 = TrialData.onset(st_TDidx);
                        t3 = TrialData.offset(st_TDidx);
                        Durations(2) = mode(diff([t2 t3],1,2));
                        t3 = t2 + Durations(2);
                        Durations(1) = mode(diff([(TrialData.onset(pst_TDidx)) t2],1,2));
                        t1 = t2 - Durations(1);
                        tstep = 100;
                        t_win = 1:tstep:1001;
                        t_FF = [-1000:10:999]';
                        
                        
                        %% PREPARE PLOT
                        
                        if PLOT_ALL
                            
                        % Make figure
                        hf = figure;
                        set(hf,'Position',figsize1)%,'NextPlot','add')
                        hold on
                        
                        pstcolors = colors(unique(TrialData(pst_TDidx,:).trID),:);
                        wincolors = flipud(winter(numel(t_win)-1));
                        
                        end %PLOT_ALL
                        
                        
                        % Preallocate
                        legstr = cell(1,numel(unique(TrialData(pst_TDidx,:).trID))); clear ip
                        vars = whos;
                        cellfun(@clear,({vars(~cellfun(@isempty,regexp({vars.name},'raster_*'))).name}))
                        stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        STIM_psth_wins = nan( numel(unique(TrialData(pst_TDidx,:).trID)), tstep, numel(t_win)-1 );
                        IR_psth_wins = [];
%                         FF       = nan(numel(unique(TrialData(pst_TDidx,:).trID)),numel(t_FF));                        
                        FF_mean  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),numel(t_FF));
                        FF_var   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),numel(t_FF));
                        r_rms_context   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        p_rms_context   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        r_drms_context  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        p_drms_context  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        
                        
                        %% GET DATA AND PLOT
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        ExcludeStim = 0;
                        Ntrials = nan(1,numel(unique(TrialData(pst_TDidx,:).trID)));
                        
                        jt=0;
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            pst_idx = find(pstid==unique(TrialData(pst_TDidx,:).trID)');
                            
                            % Make raster vectors for this transition
                            eval( sprintf('raster_x_%i = [];',pst_idx) )
                            eval( sprintf('raster_y_%i = [];',pst_idx) )
                            
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            FF_bincounts = nan(numel(trans_TDidx),numel(t_FF));
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                
                                jt=jt+1;
                                
                                psth(pst_idx,:,it) = ...
                                    Stream_FRsmooth( t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );
                                
                                stim(pst_idx,:,it) = ...
                                    SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );                                
                                
                                sp=[]; sp = spiketimes( spiketimes>=t1(trans_TDidx(it)) ...
                                    & spiketimes<t3(trans_TDidx(it)) ) - t2(trans_TDidx(it)) - 1;
                                
                                eval( sprintf('raster_x_%i = [raster_x_%i sp];',pst_idx,pst_idx) )
                                eval( sprintf('raster_y_%i = [raster_y_%i it*ones(1,numel(sp))];',pst_idx,pst_idx) )
                                
                                FF_bincounts(it,:) = sum(sp>=t_FF & sp<(t_FF+tstep),2)';
                                
                                stimTr_nSpk_minD(jt,stid) = sum(sp>=0 & sp<minDur);
                                stimTr_FR_Norm(jt,stid)   = sum(sp>=0)/(Durations(2)/1000);
                                eachTr_nSpk_minD(end+1,1) = sum(sp>=0 & sp<minDur);
                                eachTr_nSpk_Norm(end+1,1) = sum(sp>=0 & sp<Durations(2))/(Durations(2)/1000);
                                eachTr_stid_iti(end+1,1)  = stid + 0.5*round(ITIflag(is));
                                eachTr_pstid(end+1,1)     = pstid;
                                
                            end %it
                            
                            % Save N trials
                            Ntrials(1,pst_idx) = it;
                            
                            % Save the psth data for each time window
                            for iw = 2:numel(t_win)
                                IR_psth_wins(pst_idx,:,iw-1) = mean( psth( pst_idx, (t_win(iw-1):t_win(iw))+Durations(1), :) ,3,'omitnan');
                            end
                            
                            % Collect FF data
                            FF_mean(pst_idx,:) = mean(FF_bincounts,1);
                            FF_var(pst_idx,:)  = var(FF_bincounts,1);
                            
                            % Collect psth/rms correlation data
                            tmpStim = smoothFR(mean(stim(pst_idx,Durations(1)+(1:Durations(2)),:),3,'omitnan')/1000,bs_smth);
                            tmpPsth = mean(psth(pst_idx,Durations(1)+(1:Durations(2)),:),3,'omitnan');
                            
                            [r_out,p_out] = corrcoef(tmpPsth,tmpStim);
                            r_rms_context(pst_idx,1) = r_out(1,2);
                            p_rms_context(pst_idx,1) = p_out(1,2);
                            
                            [r_out,p_out] = corrcoef(tmpPsth(2:end),diff(tmpStim));
                            r_drms_context(pst_idx,1) = r_out(1,2);
                            p_drms_context(pst_idx,1) = p_out(1,2);
                            
                            
                            
                        end %pstid (prev stim id)
                        
                        
                        meanFR_all(stid) = mean(mean(mean(psth(:,Durations(1)+1:end,:),3,'omitnan')));
                        
                        
                        if PLOT_ALL
                            
                            if sum(Ntrials>=minTrs) < 2
                                ExcludeStim = 1;
                            end
                            
                            %% Create subplots
                            clear hs ip
                            hs(1)=subplot(5,1,1);   box off
                            hs(2)=subplot(5,1,2:3); box off
                            hs(3)=subplot(5,1,4:5); box off;
                            plot([0 0],[0 ymaxval],'k--')
                            hold on
                            
                        end %PLOT_ALL
                        
                        
                        % For now, plot the difference between the first 
                        % and last prev stim conditions 
                        TestDiffSpikeCount = nan(1,numel(t_win)-1);
                        if size(IR_psth_wins,1)>1
                            for iw = 2:numel(t_win)
                                % Get cumulative count of difference in n
                                % spikes
                                TestDiffSpikeCount(1,iw-1) = sum(abs( diff(IR_psth_wins([1 size(IR_psth_wins,1)],:,iw-1)./1000,1) ));
                                
                                if PLOT_ALL
                                % Plot the space between PSTHs
                                patch( -1+[t_win(iw-1):t_win(iw) t_win(iw):-1:t_win(iw-1)] ,...
                                    [IR_psth_wins(end,:,iw-1) fliplr(IR_psth_wins(1,:,iw-1))] , wincolors(iw-1,:),...
                                    'EdgeColor','none','FaceAlpha',0.65);
                                end %PLOT_ALL
                                
                            end
                        end
                        hold off
                        
                        
                        if PLOT_ALL
                        %% Plot each PSTH, raster, and stim rms trace
                        add_y = 0;
                        
                        try
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            pst_idx = find(pstid==unique(TrialData(pst_TDidx,:).trID)');
                            
                            subplot(hs(3)); hold on
                            ip(pst_idx) = plot( -Durations(1):Durations(2), ...
                                mean(psth(pst_idx,:,:),3,'omitnan') ,...
                                'Color',pstcolors(pst_idx,:),'LineWidth',4);
                            hold off
                            legstr{pst_idx} = [Info.stim_ID_key{pstid} ', n=' num2str(sum((TrialData(pst_TDidx,:).trID==pstid)))];
                            
                            subplot(hs(1)); hold on
                            plot(-Durations(1):Durations(2), mean(stim(pst_idx,:,:),3,'omitnan'),...
                                'Color',pstcolors(pst_idx,:),'LineWidth',4)
                            hold off
                            
                            if numel(eval(['raster_y_' num2str(pst_idx)]) )>0
                                subplot(hs(2)); hold on
                                plot(eval(['raster_x_' num2str(pst_idx)]),eval(['raster_y_' num2str(pst_idx)])+add_y,...
                                    '.','MarkerSize',15,'Color',pstcolors(pst_idx,:))
                                add_y = add_y + max( eval(['raster_y_' num2str(pst_idx)]) );
                                hold off
                            end
                            
                        end %pstid (prev stim id)
                        
                        catch
                            keyboard
                        end
                        
                        
                        %% Finish plot settings
                        
                        linkaxes(hs,'x')
                        xlim([-Durations(1) Durations(2)])
                        
                        subplot(hs(1)); 
                        set(gca,'xtick',[],'ytick',[])
                        
                        subplot(hs(2))
                        set(gca,'ylim',[0 add_y+1],'xtick',[])
                        
                        subplot(hs(3)); hold on
                        ylim([0 ymaxval])
                        xlabel('Time from transition (ms)')
                        ylabel('Spikes/s')
                        legend(ip,legstr,'Location','northwest')
                        hold off
                        
                        
                        % Add title
                        if ITIflag(is)
                            stimstring = sprintf('Transition to %s (ITI) stimulus',Info.stim_ID_key{stid});
                        else
                            stimstring = sprintf('Transition to %s stimulus',Info.stim_ID_key{stid});
                        end
                        suptitle(sprintf('%s   |   %s %s ch%i clu%i (%s)   |   %idB SPL, 100-%i Hz   ',...
                            stimstring, subject,session,channel,clu,unType,spl,lpn))
                                                
                        
                        end %PLOT_ALL
                        
                        
                        
                        
                        if PLOT_FF && ~ExcludeStim
                            
                            % Make figure
                            hff = figure;
                            set(hff,'Position',figsize2)
                            hold on
                            
                            pstcolors = colors(unique(TrialData(pst_TDidx,:).trID),:);
                            
                            plot([0 0],[0 5],'k--')
                            hold on
                            
                            for pstid = unique(TrialData(pst_TDidx,:).trID)'
                                
                                ip(pstid==unique(TrialData(pst_TDidx,:).trID)') = plot( t_FF+tstep/2, ...
                                    FF_var(pstid==unique(TrialData(pst_TDidx,:).trID)',:)./FF_mean(pstid==unique(TrialData(pst_TDidx,:).trID)',:) ,...
                                    'Color',pstcolors(pstid==unique(TrialData(pst_TDidx,:).trID)',:),'LineWidth',4);                                
                                legstr_ff{pstid==unique(TrialData(pst_TDidx,:).trID)'} = [Info.stim_ID_key{pstid} ', n=' num2str(sum((TrialData(pst_TDidx,:).trID==pstid)))];
                                
                                
                            end %pstid (prev stim id)
                            
                            xlim([-Durations(1) Durations(2)])
                            ylim([0 5])
                            xlabel('Time from transition (ms)')
                            ylabel('FanoFactor (100 ms bins)')
                            
                            title(sprintf('%s   |   %s %s ch%i clu%i (%s)   |   %idB SPL, 100-%i Hz   ',...
                                stimstring, subject,session,channel,clu,unType,spl,lpn))
                            
                            
                        end %PLOT_FF
                        
                        
                        
                        %++++++++++++++++++++++++++++++++++++++++++++++++++
                        %% FINISH Z-SCORE ANALYSIS
                        %++++++++++++++++++++++++++++++++++++++++++++++++++
                        
                        if CalcZ
                        % Skip stimuli that don't have at least 10 trials,
                        % or that don't have exactly 2 preceding stimulus
                        % types, for now
                        if (numel(unique(TrialData(pst_TDidx,:).trID))==2) && all(Ntrials>=minTrs)
                        
                        
                        % Bootstrap random combinations of trials
                                                
                        DistributionDiffSpikeCount = nan(nIterations,numel(t_win)-1);
                        
                        for iteration = 1:nIterations
                            
                            % Permute trials to make random assignments
                            randomtrs = pst_TDidx(randperm(length(pst_TDidx)));
                            
                            % Set up some empty vars
                            cellfun(@clear,({vars(~cellfun(@isempty,regexp({vars.name},'raster_*'))).name}))
                            stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                            psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                            IR_psth_wins = [];
                            
                            for ihalf = 1:2
                                
                                % Make raster vectors for this transition
                                eval( sprintf('raster_x_%i = [];',ihalf) )
                                eval( sprintf('raster_y_%i = [];',ihalf) )
                                
                                thesetrs = randomtrs( 1:Ntrials(1,ihalf) );
                                randomtrs( 1:Ntrials(1,ihalf) ) = [];
                                
                                % Collect FR
                                for jj = 1:numel(thesetrs)
                                    
                                    psth(ihalf,:,jj) = Stream_FRsmooth(t1(thesetrs(jj)==pst_TDidx) : t3(thesetrs(jj)==pst_TDidx) );
                                    
                                    stim(ihalf,:,jj) = SoundStream( 1, t1(thesetrs(jj)==pst_TDidx) : t3(thesetrs(jj)==pst_TDidx) );
                                    
                                    sp = []; 
                                    sp = spiketimes( spiketimes>=t1(thesetrs(jj)==pst_TDidx) & spiketimes<t3(thesetrs(jj)==pst_TDidx) )...
                                        - t2(thesetrs(jj)==pst_TDidx) - 1;
                                    eval( sprintf('raster_x_%i = [raster_x_%i sp];',ihalf,ihalf) )
                                    eval( sprintf('raster_y_%i = [raster_y_%i jj*ones(1,numel(sp))];',ihalf,ihalf) )
                                    
                                end
                                
                                for iw = 2:numel(t_win)
                                    % Save the psth data for each time window
%                                     IR_psth_wins(ihalf,:,iw-1) = mean(psth(ihalf,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
                                    IR_psth_wins(ihalf,:,iw-1) = mean( psth(ihalf,(t_win(iw-1):t_win(iw))+Durations(1), :) ,3,'omitnan');
                                end %iw
                                
                                
                            end %ihalf
                            
                            
                            % Calculate the difference between the PSTHs for
                            % each time window
                            for iw = 2:numel(t_win)
                                DistributionDiffSpikeCount(iteration,iw-1) = sum(abs( diff(IR_psth_wins(:,:,iw-1)./1000,1) ));
                            end
                            
                            
                            
                        end %iteration
                        
                        
                        %% z score
                        
                        % Calculate how far the Test spike counts are from
                        % the bootstrapped distributions
                        zScores = (TestDiffSpikeCount - mean(DistributionDiffSpikeCount,1)) ./ std(DistributionDiffSpikeCount,1);
                        
                        
                        if PLOT_ALL
                        % Add to plot
                        figure(hf)
                        subplot(hs(3)); hold on
                        for iz = 1:numel(zScores)
                            text(t_win(iz)+tstep/2,ymaxval*0.9, sprintf('%0.2f',zScores(iz)),'FontSize',14)
                        end
                        end %PLOT_ALL
                        
                        
                        end %stimulus filter for zscore analysis
                        
                        end %CalcZ
                        
                        
                        if PLOT_ALL
                        %% Save figure
                        
                        savedir = fullfile(fn.processed,'Rasters',subject,session,['ch' num2str(channel)]);
                        if ExcludeUnit || ExcludeStim
                            savedir = [savedir '/excluded'];
                        end
                        if strcmp(unType,'MU')
                            savedir = [savedir '/MU'];
                        end
                        if ~exist(savedir,'dir')
                            mkdir(savedir)
                        end
                        if ITIflag(is)
                            savename = sprintf('Trans-to-%sITI_%idB_%i_%s_%s_ch%i_%i_%s',Info.stim_ID_key{stid},spl,lpn,subject,session,channel,clu,unType);
                        else
                            savename = sprintf('Trans-to-%s_%idB_%i_%s_%s_ch%i_%i_%s',Info.stim_ID_key{stid},spl,lpn,subject,session,channel,clu,unType);
                        end
                        print_eps_kp(hf,fullfile(savedir,savename))
                        
                        end %PLOT_ALL
                        
                        
                        %**************************************************
                        %% Plot last iteration of random split, if directed
                        
                        if PLOT_RND
                            
                            % Make figure
                            hfr = figure;
                            set(hfr,'Position',figsize1,'NextPlot','add')
                            hold on
                            % Create subplots
                            clear hs ip
                            hs(1)=subplot(5,1,1);   box off
                            hs(2)=subplot(5,1,2:3); box off
                            hs(3)=subplot(5,1,4:5); box off;
                            plot([0 0],[0 ymaxval],'k--')
                            hold on
                            
                            
                            for iw = 2:numel(t_win)
                                % Plot the space between PSTHs
                                patch( -1+[t_win(iw-1):t_win(iw) t_win(iw):-1:t_win(iw-1)] ,...
                                    [IR_psth_wins(end,:,iw-1) fliplr(IR_psth_wins(1,:,iw-1))] , wincolors(iw-1,:),...
                                    'EdgeColor','none','FaceAlpha',0.65);
                            end
                            
                            for ihalf = 1:2
                                
                                subplot(hs(3)); hold on
                                ip(ihalf) = plot( -Durations(1):Durations(2), ...
                                    mean(psth(ihalf,:,:),3,'omitnan') ,...
                                    'Color',pstcolors(ihalf,:),'LineWidth',4);
                                hold off
                                %                                 legstr{pstid==unique(TrialData(pst_TDidx,:).trID)'} = [Info.stim_ID_key{pstid} ', n=' num2str(sum((TrialData(pst_TDidx,:).trID==pstid)))];
                                
                                subplot(hs(1)); hold on
                                plot(-Durations(1):Durations(2), mean(stim(ihalf,:,:),3,'omitnan'),...
                                    'Color',pstcolors(ihalf,:),'LineWidth',4)
                                hold off
                                
                                add_y=0;
                                if numel(eval(['raster_y_' num2str(ihalf)]) )>0
                                    subplot(hs(2)); hold on
                                    plot(eval(['raster_x_' num2str(ihalf)]),eval(['raster_y_' num2str(ihalf)])+add_y,...
                                        '.','MarkerSize',15,'Color',pstcolors(ihalf,:))
                                    add_y = add_y + max( eval(['raster_y_' num2str(ihalf)]) );
                                    hold off
                                end
                                
                            end %ihalf
                            
                            
                            % Finish plot settings
                            
                            linkaxes(hs,'x')
                            xlim([-Durations(1) Durations(2)])
                            
                            subplot(hs(1));
                            set(gca,'xtick',[],'ytick',[])
                            
                            subplot(hs(2))
                            set(gca,'ylim',[0 add_y+1],'xtick',[])
                            
                            subplot(hs(3)); hold on
                            ylim([0 ymaxval])
                            xlabel('Time from transition (ms)')
                            ylabel('Spikes/s')
%                             legend(ip,legstr)
                            hold off
                            
                            
                            % Add title
                            stimstring = sprintf('Transition to %s stimulus - RANDOM TRIALS',Info.stim_ID_key{stid});
                            suptitle(sprintf('%s   |   %s %s ch%i clu%i (%s)   |   %idB SPL, 100-%i Hz   ',...
                                stimstring, subject,session,channel,clu,unType,spl,lpn))
                            
                            figure(hfr); hold off
                            
                            
                        end %PLOT_RND
                        
                        
                        
                        
                        
                        
                    end %ITIflag (Trial or ITI)
                end %stid (this stim id)
                
                
                % Statistical tests to look for main effect of stimulus on
                % spikes produced
                
%                 disp('This unit has p-vals of ')
%                 p_minD = anova1(stimTr_nSpk_minD,'','off')
%                 p_norm = anova1(stimTr_FR_Norm,'','off')
%                 fprintf('----------------------------------')
%                 
%                 p_norm = anova1(stimTr_FR_Norm(~all(isnan(stimTr_FR_Norm')),:))
%                 [p_norm,~,stats] = anova1(stimTr_FR_Norm(~any(isnan(stimTr_FR_Norm')),:),'','off')
%                 multcompare(stats)
%                 
%                 [p_minD,~,stats] = anova1(stimTr_nSpk_minD(~any(isnan(stimTr_nSpk_minD')),:),'','off')
%                 multcompare(stats,'CType','bonferroni')
%                 
%                 kruskalwallis(stimTr_nSpk_minD(~any(isnan(stimTr_nSpk_minD')),:))
                
                % Would also like to test for effect of previous stimulus
                % on response
                % Currenty issue is that the dimensions are not balanced -
                % only certain stimuli are presented before each stimulus
%                 [p,aaa,stats] = anovan(eachTr_nSpk_minD,{eachTr_stid_iti eachTr_pstid},'model','linear','varnames',{'stimulus','previous_stim'});
%                 results = multcompare(stats,'Dimension',[1 2]);
%                 keyboard  


            end %lpn
        end %spl
    end %clu
end %channel
% load gong.mat;
% sound(y./2);
end %function




