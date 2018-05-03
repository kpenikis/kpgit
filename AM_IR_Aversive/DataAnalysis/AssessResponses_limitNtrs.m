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


%
% units x mean (normalized) FR across stimuli
% heatmap to visualize (+ sort to look for patterns)
% one-way repeated measures anova
%
% bootstrap segments from silent period and during AM (on spout, sample
% each stimulus, etc)
% compare FR distributions (d' value)
% similar analysis, using rcorr
%

global fn 

fn = set_paths_directories('','',1);

rng('shuffle')
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
runRCorr =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  2;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!


%%  Prepare the Resp table

if nargin<1
% Make empty data table
Resp = table;
Resp.Subject  = ' ';
Resp.Session  = ' ';
Resp.ch       = 0;
Resp.clu      = 0;
Resp.unType   = ' ';
Resp.unN      = 0;
Resp.spl      = 0;
Resp.lpn      = 0;
Resp.raw      = zeros(1,8);
Resp.nrm      = zeros(1,8);
Resp.IR_pred  = 0;
Resp.AudResp_dp  = 0;
Resp.AudResp_rc  = 0;
Resp.best_data   = cell(1,1);
Resp.bestsws     = 0;

N=0;

else  %if adding a unit to pre-existing table...
    if SUonly==1
        filename = 'Resp_allSU';
    elseif SUonly==0
        filename = 'Resp_allUnits';
    end
    Resp = readtable(fullfile(fn.processed,filename));
    
    if sum( strcmp(Resp.Subject,select_subject) &...
            strcmp(Resp.Session,select_session) &...
            Resp.ch == select_channel &...
            Resp.clu == select_clu)                   > 0
        warning('This unit already exists in the Resp table.')
        return
    end
    
    N = max(Resp.unN);
    
end

allUnResp_raw = nan(0,8);
allUnResp_nrm = nan(0,8);



%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    subjects = {select_subject};
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
    Sessions = {select_session};
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

% KS
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],24414)
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],Info.fs)


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
                
%                 fprintf(' analyzing ch %i clu %i\n',channel,clu)
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
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
                elseif  sum(Ntrials < minTrs)>1
                    keyboard
                end
                
                
                % Adjust nSnips according to min Ntrials
                if min(Ntrials)<nSnips
                    nSnips = min(Ntrials);
                    SilDATA_spk = SilDATA_spk(1:nSnips,:);
                    SilDATA_FR  = SilDATA_FR(1:nSnips,:);
                    SilDATA_z   = SilDATA_z(1:nSnips,:);
                end
                
                
%                 % Now can create cell array to save starttimes of every
%                 % trial for each stimulus (silence in column 9)
%                 %   FOR NOW: limiting # trials to minimum of all stim
%                 AllDataSnipStarttimes = cell(1,9);
%                 AllDataSnipStarttimes{1,9} = snips(idxSnips(1:nSnips))';
                
                
                
                %% Collect this unit's data
                
                dataUn_raw = nan(1,8);
                dataUn_nrm = nan(1,8);
                AMDATA_spk = nan(nSnips,minDur,8);
                AMDATA_FR  = nan(nSnips,minDur,8);
                AMDATA_z   = nan(nSnips,minDur,8);
                
                for istim = allStim'
                    
                    TDidx = all_TDidx(TrialData.trID(all_TDidx)==istim);
                    
                    % Get timestamps of onsets and offsets
                    clear t2 t3 Duration t_win
                    t2 = TrialData.onset(TDidx);
                    t3 = TrialData.offset(TDidx);
                    Duration = mode(diff([t2 t3],1,2));
                    t3 = t2 + Duration;
                    
                    % Save starttimes for rcorr analysis
                    AllDataSnipStarttimes{1,istim} = t2;
                    
                    % Collect mean FR for each trial
                    dataTr_raw = nan(nSnips,1);
                    dataTr_nrm = nan(nSnips,1);
                    nt=0;
                    its=[];
                    for it = randperm(numel(t2))
                        
                        nt=nt+1;
                        its = [its it];
                        
                        % Data collected for d' analysis
                        AMDATA_spk(nt,:,istim) = Stream_Spikes(t2(it)+[0:minDur-1]);
                        AMDATA_FR(nt,:,istim)  = Stream_FRsmooth(t2(it)+[0:minDur-1]);
                        AMDATA_z(nt,:,istim)   = Stream_zscore(t2(it)+[0:minDur-1]);
                        
                        % Data collected for comparing mean responses
                        % includes entire stimulus duration
                        dataTr_raw(nt,1) = mean(Stream_FRsmooth(t2(it):t3(it)-1));
                        dataTr_nrm(nt,1) = mean(Stream_zscore(t2(it):t3(it)-1));
                        
                        if nt==nSnips, break, end
                        
                    end %it
                    
                    % Save mean FR across trials 
                    dataUn_raw(1,istim) = mean(dataTr_raw);
                    dataUn_nrm(1,istim) = mean(dataTr_nrm);
                    
%                     % Save starttimes for rcorr analysis (if limiting ntr)
%                     AllDataSnipStarttimes{1,istim} = t2(its,1);
                    
                end %istim
                
                
                allUnResp_raw(end+1,:) = dataUn_raw;
                allUnResp_nrm(end+1,:) = dataUn_nrm;
                
                
                % -- PREDICTED IR RESPONSE --
                % Calculate prediction for IR response, based on a
                % weighted average of periodic responses
                IR_Prediction = sum( dataUn_nrm(2:6) .* ((1./AMrates)/sum(1./AMrates)) );
                
                
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                % FR d' compared to Silence
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                FRmat = format_FRmat(SilDATA_spk,AMDATA_spk);
                dprime_mat = calculate_dprime_formula(FRmat);
                
                if ~any(abs(dprime_mat(:,2))>1)
                    AudResp_dprime = 0;
%                     fprintf('EXCLUDE ch %i clu %i according to d''\n',channel,clu)
                else
                    AudResp_dprime = 1;
%                     fprintf('keep ch %i clu %i according to d''\n',channel,clu)
                end
%                 dprime_mat
                
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                % RCORR ANALYSIS
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                if runRCorr==1
                    fprintf('calculating RCORR vals for ch %i clu %i...\n',channel,clu)
                    [AudResp_rcorr,best_data,best_sw] = RCorr(Stream_Spikes,AllDataSnipStarttimes,minDur,subject,session,channel,clu);
                else
                    AudResp_rcorr=[];
                    best_data=[];
                    best_sw=[];
                end
                
                
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                % Save data to table 
                % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                Resp_addrow = { subject session channel clu unType N spl lpn dataUn_raw dataUn_nrm IR_Prediction AudResp_dprime AudResp_rcorr {best_data} best_sw };
                
                Resp = [Resp; Resp_addrow];
                
                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects


%% Finish and save Resp table

if nargin>0
    Resp(1,:) = [];
end

if SUonly==1
    filename = 'Resp_allSU';
elseif SUonly==0
    filename = 'Resp_allUnits';
end
writetable(Resp,fullfile(fn.processed,filename));



%% Play with data if wanted 

allUnResp_NRMraw = allUnResp_raw./repmat(max(allUnResp_raw,[],2),1,size(allUnResp_raw,2));
figure;
imagesc(allUnResp_NRMraw)
colormap('bone')

figure;
imagesc(allUnResp_nrm)
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
imagesc(allUnResp_raw(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('mean FR resp (Hz)\nsorted by: %s',sortBy))

hf2 = figure;
set(hf2,'Position',figsize2,'NextPlot','add')
imagesc(allUnResp_nrm(idx,2:end))
colormap('bone')
colorbar
ylabel('Unit #')
xlabel('Stimulus')
title(sprintf('norm. mean resp, (z-FR)\nsorted by: %s',sortBy))


%% Identify units to exclude

[Resp(idx,1:4) table(eval(sortBy))]
figure; hist(abs(eval(sortBy)),10)



end %function




