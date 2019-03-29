% function assess_this_unit



%% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);
if numel(dBSPL)>1 || numel(LP)>1
    keyboard
elseif isempty(dBSPL) || isempty(LP)
    warning('Not many trials in this session.')
end
AMrates = [2 4 8 16 32];



%% Label response type

fprintf('\n  %s  %s  ch %i  clu %i\n',subject,session,channel,clu)

% result = input('Unit''s response type \n 0: n/a \n 1: sparse \n 2: sustained \n 3: gap \n  ');
% switch result
%     case 0
%         RespType = 'na';
%     case 1
%         RespType = 'sparse';
%     case 2
%         RespType = 'sustained';
%     case 3
%         RespType = 'gap';
% end
        RespType = 'na';



%% Get FR over session

bs_hist = 1;
bs_smth = 20;

[Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(spiketimes,...
    length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');


% Skip if unit has very low overall firing rate
if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
    fprintf('    skipping: low FR\n')
    return
end

% Otherwise, you have a valid datapoint
N = N+1;


%% Collect data snippets from silent period

SilON   = TrialData.onset(1);
SilOFF  = TrialData.offset(1);
snipDur = 1000; %ms

% Check the duration of silence at the beginning
nSnips = floor((SilOFF - SilON)/snipDur);
if nSnips<minTrs, keyboard, end
%         fprintf('nSnips set to %i\n',nSnips)


% Randomly select snippets of the activity from the silent period
% in the beginning
snips = SilON:snipDur:SilOFF;
idxSnips = randperm(length(snips)-1);


% for rcorr, need raw spketimes
SilDATA_spk = Stream_Spikes(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
SilDATA_FR  = Stream_FRsmooth(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
SilDATA_z   = Stream_zscore(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
%         figure; hold on
%         plot(mean(SilDATA_FR,1))
%


%% Get stimulus response data

% Step through each combo of dBSPL, HP, AMdepth
for spl = dBSPL
    for lpn = LP
        
        %% Collect this unit's data
        
        % Get all stimuli presented with these parameters
        % trials without diruptive artifact while the animal was drinking
        
        [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
        if minDur~=snipDur
            keyboard
        end
        
        allStim = unique(TrialData.trID(all_TDidx));
        
        % Remove stimuli with too few trials
        if sum(Ntrials < minTrs)==1
            all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
            allStim(Ntrials<minTrs)  = [];
            Ntrials(Ntrials<minTrs) = [];
        elseif  sum(Ntrials < minTrs)>1 %if more than one stim has few trials
            keyboard
        end
        
        %                 % Adjust nSnips according to min Ntrials
        %                 if min(Ntrials)<nSnips
        %                     nSnips = min(Ntrials);
        %                     SilDATA_spk = SilDATA_spk(1:nSnips,:);
        %                     SilDATA_FR  = SilDATA_FR(1:nSnips,:);
        %                     SilDATA_z   = SilDATA_z(1:nSnips,:);
        %                 end
        %                 nSnips = min(Ntrials);
        
        
        %% First, calculate mean phase
        %  and adjust spiketimes accordingly
        try
        [IntegrationTime_spk, IntegrationTime_gap] = calculateIntegrationTime(spiketimes,TrialData,all_TDidx,Stream_Spikes,Stream_FRsmooth,AMrates,subject,session,channel,clu,RespType);
        catch
            IntegrationTime_spk = 0;
            IntegrationTime_gap = 0;
        end
        if isnan(IntegrationTime_spk)
            IntegrationTime_spk = 0;
        end
        spiketimes = spiketimes - round(IntegrationTime_spk);
        
        Stream_FRsmooth = Stream_FRsmooth(1+round(IntegrationTime_spk):end);
        Stream_zscore   = Stream_zscore(1+round(IntegrationTime_spk):end);
        Stream_Spikes   = Stream_Spikes(1+round(IntegrationTime_spk):end);
        
        %                 SpoutStream     = SpoutStream(1+round(IntegrationTime_spk):end);
        %                 SoundStream     = SoundStream(1+round(IntegrationTime_spk):end);
        
        
        
        %% Now collect rest of this unit's data
        
        Spks_per_tr = nan(max(Ntrials),1,8);
        FRtrs_raw = nan(max(Ntrials),8);
        FR_nrm    = nan(1,8);
        VSdata_spk    = nan(3,8);
        VSdata_gap    = nan(3,8);
        MeanPhase_spk = nan(1,8);
        MeanPhase_gap = nan(1,8);
        ntrs      = nan(1,8);
        
        for istim = allStim'
            
            % Skip ITI stimuli
            ITIflag = 0;
            for iiti = 1:numel(ITIflag)
                
                if istim>1 && istim<7
                    period = 1000/AMrates(istim-1);
                end
                
                TDidx = all_TDidx( TrialData.trID(all_TDidx)==istim & TrialData.ITIflag(all_TDidx) == ITIflag(iiti) );
                
                % Get timestamps of onsets and offsets
                clear t2 t3 Duration t_win
                t2 = TrialData.onset(TDidx);
                t3 = TrialData.offset(TDidx);
                Duration = mode(diff([t2 t3],1,2));
                t3 = t2 + Duration;
                
                % Preallocate
                dataTr_nrm = nan(numel(t2),1);
                raster = nan(numel(t2),Duration);
                PSTH   = nan(numel(t2),Duration);
                Gaptimes = [];
                Spktimes = [];
                
                % Collect data for each trial
                nt=0; VStr_limit = 20;
                for it = randperm(numel(t2))
                    nt=nt+1;
                    
                    % Collect spiking data (raw and normalized)
                    raster(it,:) = Stream_Spikes( (t2(it)+1): t3(it) );
                    PSTH(it,:) = Stream_FRsmooth( (t2(it)+1): t3(it) );
                    dataTr_nrm(it,1) = mean(Stream_zscore( (t2(it)+1): t3(it) ));
                    
                    % And collect spike data for VS calculation,
                    % for periodic stimuli only for now
                    if istim>1 && istim<7  %&& nt<=VStr_limit
                        sp=[]; sp = spiketimes( spiketimes>t2(it) & spiketimes<=t3(it) ) - t2(it);
                        Spktimes = [Spktimes sp];
                        % Check that 2 spike tracking methods match
                        if ~all(find(raster(it,:))==sp)
                            keyboard
                        end
                        % Collect inverted spike train
                        raster_inv = ones(size(raster(it,:)));
                        raster_inv(sp) = 0;
                        Gaptimes = [Gaptimes find(raster_inv)];
                    end
                    
                end %it
                
                ntrs(1,istim) = nt;
                
                if any(any(isnan(raster)))
                    keyboard
                end
                
                % Store data for d' FR analysis, taking mean to avoid differences in duration
                Spks_per_tr(1:size(raster,1),1,istim) = mean(raster,2);
                
                % Get mean FR across trials
                FRtrs_raw(1:size(raster,1),istim) = mean(raster,2)*1000;
                FR_nrm(1,istim) = mean(dataTr_nrm);
                
                % Get VS data
                if (istim>1 && istim<7)
                    
                    % Calculate VS and mean phase for spikes
                    [VSdata_spk(1,istim),VSdata_spk(2,istim),VSdata_spk(3,istim)] = vectorstrength(Spktimes,period);
                    MeanPhase_spk(1,istim) = meanphase(sort(Spktimes),period);
                    if MeanPhase_spk(1,istim)<0, keyboard, end
                    
                    %                         % Simulate Spktimes, to check if VS values match
                    %                         Spktimes_sim_bin = poissrnd(repmat(mean(PSTH,1),size(PSTH,1),1)/1000);
                    %                         Spktimes_sim=[];
                    %                         for it=1:size(PSTH,1)
                    %                             Spktimes_sim = [Spktimes_sim find(Spktimes_sim_bin(it,:))];
                    %                         end
                    %                         [VSdata_sim(1,istim),VSdata_sim(2,istim),VSdata_sim(3,istim)] = vectorstrength(Spktimes_sim,period);
                    %                         MeanPhase_sim(1,istim) = meanphase(sort(Spktimes_sim),period);
                    %
                    
                    % Simulate "Gaptimes", the inverted firing pattern
                    PSTH_inv = -1*mean(PSTH,1) + min(mean(PSTH,1)) + max(mean(PSTH,1));
                    Gaptimes_sim_bin = poissrnd(repmat(PSTH_inv,size(PSTH,1),1)/1000);
                    Gaptimes_sim=[];
                    for it=1:size(PSTH,1)
                        Gaptimes_sim = [Gaptimes_sim find(Gaptimes_sim_bin(it,:))];
                    end
                    
                    allUn_GR_raw(N,istim) = mean(mean(Gaptimes_sim_bin,2,'omitnan')*1000);
                    
                    % Calculate VS and mean phase for gaps (inverted
                    % but balanced firing pattern)
                    [VSdata_gap(1,istim),VSdata_gap(2,istim),VSdata_gap(3,istim)] = vectorstrength(Gaptimes_sim,period);
                    MeanPhase_gap(1,istim) = meanphase(sort(Gaptimes_sim),period);
                    if MeanPhase_gap(1,istim)<0, keyboard, end
                end
                
            end %iti
        end %istim
        
        
        % Add summary data to population vector
        allUn_FR_raw(end+1,:) = mean(FRtrs_raw,1,'omitnan');
        allUn_FR_nrm(end+1,:) = FR_nrm;
        
        
        
        %%  d' Analyses
        
        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        % FR d' compared to Silence
        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        
        % First, limit the number of trials to min of all stim
        %  (only periodic right now)
        AMDATA_spk = Spks_per_tr(:,:,~all(isnan(Spks_per_tr),1)); %2:6
        AMDATA_spk = AMDATA_spk(sum(isnan(AMDATA_spk),3)==0,:,:);
        if size(AMDATA_spk,1)<minTrs, keyboard, end
        
        FRmat = format_FRmat(SilDATA_spk,AMDATA_spk);
        dprime_mat = calculate_dprime_formula(FRmat);
        
        
        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        % Stimuli FR distribution stats
        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        % performing to be similar to Malone 2015:
        % non-parametric, and Periodic only
        
        FRdata = FRtrs_raw(:,2:6);
        FRdata = FRdata(sum(isnan(FRdata),2)==0,:,:);
        
        kw_p = kruskalwallis(FRdata,[],'off');
        
        [~,ist_max] = max(median(FRdata,1));
        [~,ist_min] = min(median(FRdata,1));
        
        wx_p = ranksum(FRdata(:,ist_min),FRdata(:,ist_max));
        
        
        
        %% Calculate integration time
        
        % For SPIKES
        
        % Wrap latency around to next cycle where necessary
        PdcPhase_spk = MeanPhase_spk(2:6);
        resets = 1+find(diff(PdcPhase_spk)<0);
        for ii = 1:numel(resets)
            PdcPhase_spk(resets(ii):end) = PdcPhase_spk(resets(ii):end) + 2*pi;
        end
        %
        %                 % Plot
        %                 figure(hf); hold on
        %                 subplot(hs_spk); hold on
        %                 plot(AMrates,rad2deg(PdcPhase_spk),'-','LineWidth',1.5)
        %                 hold off
        %
        %                 % Fit slope to find integration time
        %                 weights_spk = allUn_FR_raw(end,2:6) .* VSdata_spk(1,2:6);
        %
        %                 [IntegrationTime_nt,~,IT_mse] = lscov(AMrates',rad2deg(PdcPhase_spk'),ntrs(2:6));
        %                 [IntegrationTime_spk,~,IT_mse] = lscov(AMrates',rad2deg(PdcPhase_spk'),weights_spk);
        %
        %                 if any(rad2deg(PdcPhase_spk)<0)
        %                     keyboard
        %                 end
        
        
        % For GAPS
        
        % Wrap latency around to next cycle where necessary
        PdcPhase_gap = MeanPhase_gap(2:6);
        resets = 1+find(diff(PdcPhase_gap)<0);
        for ii = 1:numel(resets)
            PdcPhase_gap(resets(ii):end) = PdcPhase_gap(resets(ii):end) + 2*pi;
        end
        
        %                 % Plot
        %                 subplot(hs_gap); hold on
        %                 plot(AMrates,rad2deg(PdcPhase_gap),'-','LineWidth',1.5)
        %                 hold off
        %
        %                 subplot(hs_diff); hold on
        %                 plot(AMrates,rad2deg(PdcPhase_gap)-rad2deg(PdcPhase_spk),'-','LineWidth',1.5)
        %                 hold off
        %
        %                 % Fit slope to find integration time
        %                 weights_gap = allUn_GR_raw(end,2:6) .* VSdata_gap(1,2:6);
        %
        %                 [IntegrationTime_gap,~,IT_mse] = lscov(AMrates',rad2deg(PdcPhase_gap'),weights_gap);
        %
        %                 if any(rad2deg(PdcPhase_gap)<0)
        %                     keyboard
        %                 end
        
        
        %                 % Plot response phases for this unit
        %
        %                 hf_tmp=figure; hold on
        %                 ip(1)=plot(AMrates,rad2deg(PdcPhase_spk),'-k','LineWidth',2.5);
        %                 ip(2)=plot(AMrates,rad2deg(PdcPhase_gap),'-r','LineWidth',2.5);
        %                 for ir = 1:5
        %                     plot(AMrates(ir),rad2deg(PdcPhase_spk(ir)),'.k','MarkerSize',10*weights_spk(ir))
        %                     plot(AMrates(ir),rad2deg(PdcPhase_gap(ir)),'.r','MarkerSize',10*weights_gap(ir))
        %                 end
        %                 ylim([0 max( rad2deg(max(PdcPhase_spk)+pi/2), rad2deg(max(PdcPhase_gap)+pi/2)) ])
        %                 xlim([0 34])
        %                 axis square
        % %                 title(sprintf('Integration Time \nspks: %0.1f ms, \\color{red}gaps: %0.1f ms',IntegrationTime_spk,IntegrationTime_gap)...
        % %                     ,'Interpreter','tex')
        %                 title(sprintf('%s %s ch%i clu%i',subject,session,channel,clu),'Interpreter','none')
        %                 legend(ip,{sprintf('spks: %0.1f ms',IntegrationTime_spk) sprintf('gaps: %0.1f ms',IntegrationTime_gap)},'Location','southeast')
        %                 set(gca,'xtick',AMrates)
        %                 xlabel('AM rate (Hz)')
        %                 ylabel('Unwrapped mean phase (deg)')
        %                 hold off
        %
        
        
        %% Save data to structure
        
        UnitData(N).Subject     = subject;
        UnitData(N).Session     = session;
        UnitData(N).Channel     = channel;
        UnitData(N).Clu         = clu;
%         UnitData(N).unType      = unType;
        UnitData(N).spl         = spl;
        UnitData(N).lpn         = lpn;
        UnitData(N).BaseFR      = sum(Stream_Spikes(SilON:SilOFF))/(SilOFF-SilON)*1000;
        UnitData(N).FR_raw_tr   = FRtrs_raw;
        UnitData(N).kw_p        = kw_p;
        UnitData(N).wx_p        = wx_p;
        UnitData(N).VSdata_spk  = VSdata_spk;
        UnitData(N).VSdata_gap  = VSdata_gap;
        UnitData(N).Phase_spk   = rad2deg(PdcPhase_spk);
        UnitData(N).Phase_gap   = rad2deg(PdcPhase_gap);
        UnitData(N).IntTime_spk = IntegrationTime_spk;
        UnitData(N).IntTime_gap = IntegrationTime_gap;
        UnitData(N).ntr         = ntrs(2:6);
        UnitData(N).FR_nrm      = FR_nrm;
        UnitData(N).dp_mat      = dprime_mat;
        
        
        %% Save info to table
        
        add_row = { subject session channel clu RespType };
        
        UnitInfo = [ UnitInfo; add_row ];
        
        %
        %                 % Save fig
        % %                 savedir = fullfile(fn.processed,'Units',sprintf('%s_%s_%i_%i',subject,session,channel,clu));
        %                 savedir = fullfile(fn.processed,'Units',RespType);
        %                 if ~exist(savedir,'dir')
        %                     mkdir(fullfile(savedir,'eps'))
        %                     mkdir(fullfile(savedir,'svg'))
        %                 end
        %                 print_eps_kp(hf_tmp,fullfile(savedir,'eps',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))
        %                 print_svg_kp(hf_tmp,fullfile(savedir,'svg',sprintf('RespPhase_%s_%s_%i_%i',subject,session,channel,clu)))
        
        % Close this unit's fig
        %                 close(hf_tmp)
        
    end %lpn
end %spl


