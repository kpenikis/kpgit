% function assess_this_unit


%% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);
% if numel(dBSPL)>1 || numel(LP)>1
%     keyboard    %if WWWlf_395 DA: dBSPL = 75;
% end


%% Label response type

fprintf('  %s  %s  ch %i  clu %i\n',subject,session,channel,clu)

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

[Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(unique(round(spiketimes)),...
    length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');


% Skip if unit has very low overall firing rate
if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
    fprintf('    skipping: low FR\n')
    return
end


%% Collect data snippets from silent period

SilON   = TrialData.onset(1);
SilOFF  = TrialData.offset(1);
snipDur = 1000; %ms

% Check the duration of silence at the beginning
nSnips = floor((SilOFF - SilON)/snipDur);
% if nSnips<15, keyboard, end


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
for spl = dBSPL'
    for lpn = LP'
        
        %% Collect this unit's data
        
        % Get all stimuli presented with these parameters
        % trials without diruptive artifact while the animal was drinking
        
        [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn,1);
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
            continue
        end
        
        %                 % Adjust nSnips according to min Ntrials
        %                 if min(Ntrials)<nSnips
        %                     nSnips = min(Ntrials);
        %                     SilDATA_spk = SilDATA_spk(1:nSnips,:);
        %                     SilDATA_FR  = SilDATA_FR(1:nSnips,:);
        %                     SilDATA_z   = SilDATA_z(1:nSnips,:);
        %                 end
        %                 nSnips = min(Ntrials);
        
        
        
        % Otherwise, you have a valid datapoint
        N = N+1;
        
        
        %% First, calculate Integration Time
        %  and adjust spiketimes accordingly
        
        try
        IntegrationTime_spk = calculateIntegrationTime(unique(round(spiketimes)),TrialData,all_TDidx,Stream_Spikes,Stream_FRsmooth,AMrates,subject,session,channel,clu,RespType);
        catch
            warning('integration time program didnt run properly')
            IntegrationTime_spk = 0;
        end
        if isnan(IntegrationTime_spk)
            IntegrationTime_spk = 0;
        end
        spiketimes = unique(round(spiketimes - IntegrationTime_spk));
        
        % Re-run with shifted spiketimes
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');
        
        
        %% Now collect rest of this unit's data
        %   AFTER shifting spike times!
        
        Spks_per_tr = nan(max(Ntrials),1,8);
        FRtrs_raw = nan(max(Ntrials),8);
        FR_nrm    = nan(1,8);
        VSdata_spk    = nan(3,8);
        MeanPhase_spk = nan(1,8);
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
        
        
        
        %% Save data to structure
        
        UnitData(N).Subject     = subject;
        UnitData(N).Session     = session;
        UnitData(N).Shank       = shank;
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
        UnitData(N).Phase_spk   = rad2deg(MeanPhase_spk);
        UnitData(N).IntTime_spk = IntegrationTime_spk;
        UnitData(N).ntr         = ntrs(2:6);
        UnitData(N).FR_nrm      = FR_nrm;
        UnitData(N).dp_mat      = dprime_mat;
        
        
        %% Save info to table
        
        add_row = { subject session channel clu RespType };
        [add_row{end+1:size(UnitInfo,2)}] = deal(nan);
        
        UnitInfo = [ UnitInfo; add_row ];
        
        
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


