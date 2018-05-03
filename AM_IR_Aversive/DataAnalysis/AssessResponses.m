function AssessResponses(select_subject, select_session, select_channel, select_clu)
%
%  AssessResponses( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   Uses the TrialData (newer) version of saving stimulus info.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%
%   If no input arguments, runs through all units, overwriting previous
%   file. If a unit is specified in the input arguments, program adds the
%   unit to the existing table. If the unit entered as the input already
%   exists in the table, the program will exit. In the future, could allow
%   the row to update individually.
%
%  KP, 2018-03
%


global fn 

fn = set_paths_directories('','',1);

rng('shuffle')
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
runRCorr =  1;  nIterations = 1000;
%!!!!!!!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!


%%  Prepare the Resp table

if nargin<1 %start over, re-running all units
    
    RERUN = input('If want to re-run ALL UNITS, enter 1. Otherwise, enter 0.   ');
    
end

if RERUN == 1
    
    Resp = struct;
    N=0;

else  %to add units to pre-existing struct
    
    if SUonly==1
        filename = 'RespStruct_allSU';
    elseif SUonly==0
        filename = 'RespStruct_allUnits';
    end
    
    q = load(fullfile(fn.processed,filename));
    Resp = q.Resp;
    
%     if sum( strcmp(Resp.Subject,select_subject) &...
%             strcmp(Resp.Session,select_session) &...
%             Resp.ch == select_channel &...
%             Resp.clu == select_clu)                   > 0
%         warning('This unit already exists in the Resp table.')
%         return
%     end
    
    N = numel(Resp);
    
end

allUn_FR_raw = nan(0,8);
allUn_FR_nrm = nan(0,8);



%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    if ~iscell(select_subject)
        subjects = {select_subject};
    end
else
    subjects = {'WWWf_253400' 'WWWlf_253395'};
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    
    switch subject
        case 'WWWf_253400'
            subjcol = 'k';
        case 'WWWlf_253395'
            subjcol = 'k';
    end
     

%%  SESSIONS

% Get list of sessions to check for sorted data

if nargin>1 && exist('select_session','var')
    if ~iscell(select_session)
        Sessions = {select_session};
    end
else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
    Sessions = [];
    for ifn = 1:numel(SpkFns)
        if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
            Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
        end
    end
    
end
% Sessions = flipud(Sessions);

% Step through each session
for sess = Sessions'
    
session = char(sess);

%%
% Load data files
fn = set_paths_directories(subject,session,1);
fprintf('Loading sess %s...\n',session)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));


%% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);
AMrates = [2 4 8 16 32];


%% STEP THROUGH EACH CHANNEL

if nargin>2 && exist('select_channel','var')
    Channels = select_channel;
else
    Channels =  [1:7 9:16];
end

for channel = Channels
        
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3) %no valid clus for this channel
            continue
        else
            Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
        end
    else
        Clus = select_clu;
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = Clus'
        
        
        % Skip if clu already in struct, unless you want to rerun all
        if RERUN == 0
            if sum( strcmp({Resp.Subject},subject) & strcmp({Resp.Session},session) & [Resp.Channel]==channel & [Resp.Clu]==clu ) == 1
                continue
            end
        end
        
        
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
            continue
        end
        
        
%         fprintf('-----------------------------\n')
        
        
        
        
        
        %% ################     NOW THE FUN STARTS     ##################
        
        
        
        
        
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Collect FR over session
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        bs_hist = 1;
        bs_smth = 20;
        
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');
        
        
        % Check if unit has very low overall firing rate
        if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
            continue
        end
        
        N = N+1;
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Collect data snippets from silent period
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
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
        
        % Create cell array to save starttimes, of every trial
        % for each stimulus (silence in column 9), for RCorr
        %   (not limiting # trials)
        AllDataSnipStarttimes = cell(1,9);
        AllDataSnipStarttimes{1,9} = snips(1:end-1)';
        
        %%%%%% for rcorr, need raw spketimes
        SilDATA_spk = Stream_Spikes(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
        SilDATA_FR  = Stream_FRsmooth(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
        SilDATA_z   = Stream_zscore(repmat(snips(idxSnips(1:nSnips))',1,snipDur) + repmat([0:snipDur-1],nSnips,1));
%         figure; hold on
%         plot(mean(SilDATA_FR,1))
%                 


        %%
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
                %% Collect this unit's data
                
                fprintf(' analyzing ch %i clu %i\n',channel,clu)
                
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
                
% %                 % Adjust nSnips according to min Ntrials
% %                 if min(Ntrials)<nSnips
% %                     nSnips = min(Ntrials);
% %                     SilDATA_spk = SilDATA_spk(1:nSnips,:);
% %                     SilDATA_FR  = SilDATA_FR(1:nSnips,:);
% %                     SilDATA_z   = SilDATA_z(1:nSnips,:);
% %                 end
                
%                 nSnips = min(Ntrials);
                Spks_per_tr = nan(max(Ntrials),1,8);
                FRtrs_raw = nan(max(Ntrials),8);
                FR_nrm = nan(1,8);
                VSdata = nan(3,8);
                
                
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
                    
                    % Save starttimes for rcorr analysis
                    AllDataSnipStarttimes{1,istim} = t2;
                    
                    % Preallocate
                    dataTr_nrm = nan(numel(t2),1);
                    raster = nan(numel(t2),Duration);
                    Spktimes = [];
                    
                    % Collect data for each trial
                    nt=0; VStr_limit = 20;
                    for it = randperm(numel(t2))
                        nt=nt+1;
                        
                        % Collect spiking data (raw and normalized)
                        raster(it,:) = Stream_Spikes( (t2(it)+1): t3(it) );
                        dataTr_nrm(it,1) = mean(Stream_zscore( (t2(it)+1): t3(it) ));
                                                                        
                        % And collect spike data for VS calculation,
                        % for periodic stimuli only for now
                        if istim>1 && istim<7 && nt<=VStr_limit
                            sp=[]; sp = spiketimes( spiketimes>t2(it) & spiketimes<=(t3(it)-1) ) - t2(it);
                            Spktimes = [Spktimes sp];
                        end
                        
                    end %it
                    
                    if any(any(isnan(raster)))
                        keyboard
                    end
                    
                    % Save data for d' FR analysis, taking mean to avoid differences in duration
                    Spks_per_tr(1:size(raster,1),1,istim) = mean(raster,2);
                    
                    % Save mean FR across trials 
                    FRtrs_raw(1:size(raster,1),istim) = mean(raster,2)*1000;
                    FR_nrm(1,istim) = mean(dataTr_nrm);
                    
                    % Save VS data
                    if (istim>1 && istim<7)
                        [VSdata(1,istim),VSdata(2,istim),VSdata(3,istim)] = vectorstrength(Spktimes,period);
                    end
                    
                    end %iti
                end %istim
                
                allUn_FR_raw(end+1,:) = mean(FRtrs_raw,1,'omitnan');
                allUn_FR_nrm(end+1,:) = FR_nrm;
                
                
                
                %%  Analyses to check for significnat responses/tuning
                
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
                
                
                
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                % RCORR ANALYSIS
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                if runRCorr==1
                    fprintf('calculating RCORR vals for ch %i clu %i...\n',channel,clu)
                    [best_data,best_PC,best_sw] = RCorr( ...
                        Stream_Spikes,AllDataSnipStarttimes,minDur,nIterations,...
                        subject,session,channel,clu);
                else
                    best_data=[];
                    best_PC  =[];
                    best_sw  =[];
                end
                
                
                
                
                %% Save data to structure
                
                Resp(N).Subject    = subject;
                Resp(N).Session    = session;
                Resp(N).Channel    = channel;
                Resp(N).Clu        = clu;
                Resp(N).unType     = unType;
                Resp(N).spl        = spl;
                Resp(N).lpn        = lpn;
                Resp(N).FR_raw_tr  = FRtrs_raw;
                Resp(N).kw_p       = kw_p;
                Resp(N).wx_p       = wx_p;
                Resp(N).VSdata     = VSdata;
                Resp(N).FR_nrm     = FR_nrm;
                Resp(N).dp_mat     = dprime_mat;
                Resp(N).RCorr.data = best_data;
                Resp(N).RCorr.PC   = best_PC;
                Resp(N).RCorr.sw   = best_sw;
                
                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects


%% Finish and save Resp structure

if SUonly==1
    filename = 'RespStruct_allSU';
elseif SUonly==0
    filename = 'RespStruct_allUnits';
end

save(fullfile(fn.processed,filename),'Resp','-v7.3');

return


%% Play with data if wanted 

allUnResp_NRMraw = allUn_FR_raw./repmat(max(allUn_FR_raw,[],2),1,size(allUn_FR_raw,2));
figure;
imagesc(allUnResp_NRMraw)
colormap('bone')

figure;
imagesc(allUn_FR_nrm)
colormap('bone')


% Sort according to some feature

sortBy = 'pdc_raw_resp_range';
switch sortBy
    case 'warn_raw'
        [warn_raw,idx]=sort(Resp.raw(:,1),1);
    case 'raw_32'
        [raw_32,idx]=sort(Resp.raw(:,6),1);
    case 'nrm_32'
        [nrm_32,idx]=sort(Resp.nrm(:,6),1);
    case 'IR_pred'
        [IR_pred,idx]=sort(Resp.IR_pred,1);
    case 'diff_nrm_2_32'
        [diff_nrm_2_32,idx]=sort(Resp.nrm(:,2)-Resp.nrm(:,6),1);
    case 'pdc_raw_resp_range'
        [pdc_raw_resp_range,idx]=sort(range(Resp.raw(:,2:6),2),1);
    case 'pdc_nrm_resp_range'
        [pdc_nrm_resp_range,idx]=sort(range(Resp.nrm(:,2:6),2),1);        
end


%% Plots
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4) scrsz(3)/8 scrsz(4)];
figsize2 = [scrsz(3)/8 scrsz(4) scrsz(3)/8 scrsz(4)];

hf1 = figure;
set(hf1,'Position',figsize1,'NextPlot','add')
imagesc(allUn_FR_raw(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('mean FR resp (Hz)\nsorted by: %s',sortBy))

hf2 = figure;
set(hf2,'Position',figsize2,'NextPlot','add')
imagesc(allUn_FR_nrm(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('norm. mean resp, (z-FR)\nsorted by: %s',sortBy))


%% Identify units to exclude

[Resp(idx,1:4) table(eval(sortBy))]
figure; hist(abs(eval(sortBy)),10)



end %function




