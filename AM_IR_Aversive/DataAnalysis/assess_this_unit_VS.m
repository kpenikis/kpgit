% function assess_this_unit

fprintf('  %s  %s  ch %i  clu %i\n',subject,session,channel,clu)


%% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);
if numel(dBSPL)>1 || numel(LP)>1
    keyboard   
end


%% Get FR over session

bs_hist = 1;
bs_smth = 20;

[Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(round(spiketimes),...
    length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');


% Skip if unit has very low overall firing rate
if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
    fprintf('    skipping: low FR\n')
    return
end


%% Collect data snippets from silent period

SilON   = TrialData.onset(1);
SilOFF  = TrialData.offset(1);
snipDur = 1500; %ms

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
        [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,[],spl,lpn,0);
        
        allStim = unique(TrialData.trID(all_TDidx));
        
        % Remove stimuli with too few trials
        if sum(Ntrials < minTrs)==1
            all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
            allStim(Ntrials<minTrs)  = [];
            Ntrials(Ntrials<minTrs) = [];
        elseif  sum(Ntrials < minTrs)>1 %if more than one stim has few trials
            continue
        end
        
        % Otherwise, you have a valid datapoint
        N = N+1;
        
        
        %% Now collect rest of this unit's data
        
        nStim       =  numel(allStim);
        FRtrs_raw   =  nan(max(Ntrials),nStim);
        FR_nrm      =  nan(1,nStim);
        ntrs        =  nan(1,nStim);
        Dur         =  nan(1,nStim);
        
        for istim = allStim'
            
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==istim );
            
            % Get timestamps of onsets and offsets
            clear t2 t3 Duration t_win
            t2 = TrialData.onset(TDidx);
            t3 = TrialData.offset(TDidx);
            Duration = mode(diff([t2 t3],1,2));
            t3 = t2 + Duration;
            
            % Preallocate
            dataTr_nrm  =  nan(numel(t2),1);
            raster      =  nan(numel(t2),Duration);
            PSTH        =  nan(numel(t2),Duration);
            
            % Collect data for each trial
            for it = 1:numel(t2)
                
                % Collect spiking data (raw and normalized)
                raster(it,:)     =  Stream_Spikes(      (t2(it)+1) : t3(it) );
                PSTH(it,:)       =  Stream_FRsmooth(    (t2(it)+1) : t3(it) );
                dataTr_nrm(it,1) =  mean(Stream_zscore( (t2(it)+1) : t3(it) ));
                
            end %it
            
            ntrs(1,istim) = numel(t2);
            Dur(1,istim)  = Duration;
            
            if any(any(isnan(raster)))
                keyboard
            end
            
            % Get mean FR across trials
            FRtrs_raw(1:size(raster,1),istim) = mean(raster,2)*1000;
            FR_nrm(1,istim) = mean(dataTr_nrm);
            
        end %istim
        
        
        %% Save data to structure
        
        UnitData(N).Subject     = subject;
        UnitData(N).Session     = session;
        UnitData(N).Shank       = shank;
        UnitData(N).Channel     = channel;
        UnitData(N).Clu         = clu;
        UnitData(N).spl         = spl;
        UnitData(N).lpn         = lpn;
        UnitData(N).BaseFR      = sum(Stream_Spikes(SilON:SilOFF))/(SilOFF-SilON)*1000;
        UnitData(N).FR_raw_tr   = FRtrs_raw;
        UnitData(N).ntr         = ntrs;
        UnitData(N).Dur         = Dur;
        UnitData(N).FR_nrm      = FR_nrm;
        
        
        %% Save info to table
        
        add_row = { subject session channel clu };
        UnitInfo = [ UnitInfo; add_row ];
        
    end %lpn
end %spl


