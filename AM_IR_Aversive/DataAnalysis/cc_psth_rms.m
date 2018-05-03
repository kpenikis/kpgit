function cc_psth_rms(select_subject, select_session, select_channel, select_clu)
%
%  corr_psth_rms( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   Uses the TrialData (newer) version of saving stimulus info.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%
%  KP, 2018-04
%



global ms_pad fn
ms_pad = 500;

rng('shuffle')
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  2;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
PLOT_FF = 0;

%% Set parameters 

tstep = 50;
t_FF = [-1000:10:999]';

foo = nan(0,4);
corr_ContextData = cell(8,8);
[corr_ContextData{:}] = deal(foo);


%%
% Make empty data table
Resp = table;
Resp.Subject  = ' ';
Resp.Session  = ' ';
Resp.ch       = 0;
Resp.clu      = 0;
Resp.unType   = ' ';
Resp.unN      = 0;
Resp.stim     = 0;
Resp.iti      = 0;
Resp.spl      = 0;
Resp.lpn      = 0;
Resp.r_RMS    = 0;
Resp.p_RMS    = 0;
Resp.r_dRMS   = 0;
Resp.p_dRMS   = 0;

N=0;



%% 
% Set up figure options

colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [colors; 0.7.*bone(2)];

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
figsize1 = [1 2*scrsz(4)/3 scrsz(3) 2*scrsz(4)/3];
figsize2 = [1 1 scrsz(3) scrsz(4)/3];



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

fn = set_paths_directories(subject,'',1);

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
if numel(dBSPL)>1 || numel(LP)>1, keyboard, end
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
                
                
                
                %% Collect this unit's data
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                for istim = allStim'
                    
                    %% FIRST, COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
                    
                    %%%  plot for Trial and ITI stimuli separately
                    ITIflag = unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for iti = 1%:numel(ITIflag)
                        
                        st_TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(iti));
                        
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
                        
                        
                        %% GET DATA AND PLOT
                        
                        % Preallocate
                        stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        r_rms_context   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        p_rms_context   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        r_drms_context  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        p_drms_context  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),1);
                        CompositeStim   = nan(0,Durations(2));
                        CompositePsth   = nan(0,Durations(2));
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        ExcludeStim = 0;
                        Ntrials = nan(1,numel(unique(TrialData(pst_TDidx,:).trID)));
                        
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            FF_bincounts = nan(numel(trans_TDidx),numel(t_FF));
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                
                                pst_idx = find(pstid==unique(TrialData(pst_TDidx,:).trID)');
                                
                                psth(pst_idx,:,it) = ...
                                    Stream_FRsmooth( t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );
                                
                                stim(pst_idx,:,it) = ...
                                    SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );      
                                
                                CompositeStim = [CompositeStim; SoundStream(1, t2(trans_TDidx(it)) : t3(trans_TDidx(it))-1 ) ];
                                CompositePsth = [CompositePsth; Stream_FRsmooth( t2(trans_TDidx(it)) : t3(trans_TDidx(it))-1 ) ];
                                
%                                 sp=[]; sp = spiketimes( spiketimes>=t1(trans_TDidx(it)) ...
%                                     & spiketimes<t3(trans_TDidx(it)) ) - t2(trans_TDidx(it)) - 1;
                                
%                                 FF_bincounts(it,:) = sum(sp>=t_FF & sp<(t_FF+tstep),2)';
                                
                            end %it
                            
                            Ntrials(1,pst_idx) = it;
                            
                            % Collect psth/rms correlation data, by context
                            tmpStim = smoothFR(mean(stim(pst_idx,Durations(1)+(1:Durations(2)),:),3,'omitnan')/1000,bs_smth);
                            tmpPsth = mean(psth(pst_idx,Durations(1)+(1:Durations(2)),:),3,'omitnan');
                            
                            [r_out,p_out] = corrcoef(tmpPsth,tmpStim);
                            r_rms_context(pst_idx,1) = r_out(1,2);
                            p_rms_context(pst_idx,1) = p_out(1,2);
                            
                            [r_out,p_out] = corrcoef(tmpPsth(2:end),diff(tmpStim));
                            r_drms_context(pst_idx,1) = r_out(1,2);
                            p_drms_context(pst_idx,1) = p_out(1,2);
                            
                            if it>=minTrs
                                corr_ContextData{istim,pstid}(end+1,:) = [r_rms_context(pst_idx,1) p_rms_context(pst_idx,1) r_drms_context(pst_idx,1) p_drms_context(pst_idx,1)];
                            end
                             
                        end %pstid (prev stim id)
                        
                        
                        % Find a better way to exclude individual stimuli
                        if sum(Ntrials)<minTrs
                            ExcludeStim = 1;
                        else
                            
                            [r_rms, p_rms ] = corrcoef(mean(CompositePsth,1,'omitnan'), mean(CompositeStim,1,'omitnan'));
                            [r_drms,p_drms] = corrcoef(mean(CompositePsth(:,2:end),1,'omitnan'), diff(mean(CompositeStim,1,'omitnan')));
                        
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        % Save data to table
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        Resp_addrow = { subject session channel clu unType N istim iti-1 spl lpn r_rms(1,2) p_rms(1,2) r_drms(1,2) p_drms(1,2) };
                        
                        Resp = [Resp; Resp_addrow];
                        
                        end
                        
                    end %ITIflag (Trial or ITI)                    
                end %istim
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects

corr_ContextData
Resp(1,:) = [];

figure;
hist(Resp(Resp.p_RMS<0.01&Resp.p_dRMS<0.01,:).stim,0:9)

keyboard


%% Plot population average data


% close all
clear hf hs hfb

hfb=figure;
twin = [0 500];
t_vec = t_FF'+tstep/2;
title(sprintf('Average FF for each stimulus during first %i ms',twin(2)))

for istim = 1:size(corr_ContextData,1)
    
    hf(istim) = figure;
    set(hf(istim),'Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)],'NextPlot','add')
    
    nsp = 1+sum(~cellfun(@isempty,corr_ContextData(istim,:)));
    for isp=1:nsp
        hs(isp)=subplot(nsp,1,isp);
        plot([0 0],[0.4 1.6],':k','LineWidth',0.75)
        xlim([t_FF(1) t_FF(end)]+tstep/2)
        ylim([0 3])
%         ylim([0.4 1.6])
    end
    
    stimdata = [];
    isp=0;
    for ipst = 1:size(corr_ContextData,2)
        if ~isempty(corr_ContextData{istim,ipst})
            
            isp=isp+1;
            data = corr_ContextData{istim,ipst};
            stimdata = [stimdata; data];
            
            subplot(hs(isp))
            hold on
%             plot(t_vec, mean(data,1,'omitnan')+std(data,1,'omitnan')./sqrt(size(data,1)),...
%                 'LineWidth',1,'Color',colors(ipst,:))
%             plot(t_vec, mean(data,1,'omitnan')-std(data,1,'omitnan')./sqrt(size(data,1)),...
%                 'LineWidth',1,'Color',colors(ipst,:))
%             plot(t_vec, mean(data,1,'omitnan'),'LineWidth',3,'Color',colors(ipst,:))
            plot(repmat(t_vec,size(data,1),1)', data','LineWidth',1,'Color',colors(ipst,:))
            hold off
            
            subplot(hs(end))
            hold on
            plot(t_vec, mean(data,1,'omitnan')+std(data,1,'omitnan')./sqrt(size(data,1)),...
                'LineWidth',1,'Color',colors(ipst,:))
            plot(t_vec, mean(data,1,'omitnan')-std(data,1,'omitnan')./sqrt(size(data,1)),...
                'LineWidth',1,'Color',colors(ipst,:))
            plot(t_vec, mean(data,1,'omitnan'),'LineWidth',3,'Color',colors(ipst,:))
            ylim([0.4 1.6])
            hold off
        end
    end
    
    xlabel('Time from transition (ms)')
    ylabel('FanoFactor (100 ms bins)')
    suptitle(sprintf('Transitions to %s stimulus',Info.stim_ID_key{istim}))
    hold off
    
    
    % Add to bar plot of average FF for each stim
    figure(hfb); hold on
    ib=bar(istim,mean(mean(stimdata(:,t_vec>twin(1)&t_vec<twin(2)),2,'omitnan'),'omitnan'),...
        'EdgeColor','none','FaceColor',colors(istim,:));
    if istim==1
        ib.EdgeColor = 'k';
    end
    plot([istim istim],mean(mean(stimdata(:,t_vec>twin(1)&t_vec<twin(2)),2,'omitnan'),'omitnan')+...
        [-std(mean(stimdata(:,t_vec>twin(1)&t_vec<twin(2)),2,'omitnan'),'omitnan') std(mean(stimdata(:,t_vec>twin(1)&t_vec<twin(2)),2,'omitnan'),'omitnan')]...
        ./sqrt(size(stimdata,1)),'Color','k','LineWidth',2);
    

end %istim

%%
keyboard

end %function




