function ContextMTFs(select_subject, select_session, select_channel, select_clu)
%
%  ContextMTFs( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   Uses the TrialData (newer) version of saving stimulus info.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%   Option to plot only robust responders.
%
%  KP, 2018-04
%



global fn 
close  all

%!!!!!!!!!!!!!!!!!
PLOT_EACH_MTF = 1;
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
RobustUn =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  2;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  8;
%!!!!!!!!!!!!!!!!!
spktimeSHIFT = 0;
%!!!!!!!!!!!!!!!!!

% Set up figure options
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
figsize = [1 1 scrsz(3) scrsz(4)/2];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];



%% Load Resp data table (with RCorr results)

fn = set_paths_directories('','',1);
Resp = readtable(fullfile(fn.processed,'Resp_allSU'));

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

N=0;

foo = nan(0,1);
BinarySpikeData = cell(8,8);
[BinarySpikeData{:}] = deal(foo);

% FR_GroupData=struct;
% FR_GroupData.mean = fii;
% FR_GroupData.std  = fii;

MTF_ts = [100 250 500 nan]';
foo = nan(8,8,3);

FR_Population=struct;
for ii=1:numel(MTF_ts)
    if isnan(MTF_ts(ii))
        MTF_ts_str{ii} = 'durFull';
    else
        MTF_ts_str{ii} = sprintf('dur%i',MTF_ts(ii));
    end
    FR_Population.(MTF_ts_str{ii}) = foo;
end


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


if RobustUn %skip sessions with no Robust SUs
    SessResp = Resp(strcmp(Resp.Subject,subject)&strcmp(Resp.Session,session),:);
    SessResp = SessResp(SessResp.AudResp_rc==1,:);
    if isempty(SessResp),  continue,   end
    Channels = unique(SessResp.ch)';
else 
    Channels =  [1:7 9:16];
end


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
if numel(dBSPL)>1 || numel(LP)>1
    keyboard
end
AMrates = [2 4 8 16 32];


%% STEP THROUGH EACH CHANNEL

for channel = Channels
        
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if RobustUn==0
        if nargin<4
            if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3) %no valid clus for this channel
                continue
            else
                Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
            end
        else
            Clus = select_clu;
        end
    else
        Clus = SessResp(SessResp.ch==channel,:).clu';
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
        spiketimes = spiketimes+spktimeSHIFT;
        
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
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
%                 fprintf(' analyzing ch %i clu %i\n',channel,clu)
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
                all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                
                allStim = unique(TrialData.trID(all_TDidx));
                
%                 % Remove stimuli with too few trials
%                 if sum(Ntrials < minTrs)==1
%                     all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
%                     allStim(Ntrials<minTrs)  = [];
%                     Ntrials(Ntrials<minTrs) = [];
%                 elseif  sum(Ntrials < minTrs)>1
%                     keyboard
%                 end
                
                
                %==========================================================
                % Set up empty data matrices for this unit
                %==========================================================
                for ii = 1:numel(MTF_ts)
                    eval(sprintf('%s_data = foo;',MTF_ts_str{ii}))
                end
                %==========================================================
                
                
                %% Collect this unit's data
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                
                for istim = allStim'
                    
                    %% 1) COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
                    
                    %%%  plot for Trial and ITI stimuli separately
                    ITIflag = unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for is = 1:numel(ITIflag)
                        
                        st_TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
                        
                        pst_TDidx = nan(size(st_TDidx));
                        skip_it = [];
                        
                        for it = 1:numel(st_TDidx)
                            
                            % Get TD index of previous (clean) trial
                            
                            if (find(all_TDidx==st_TDidx(it))-1)==0 %if very first trial
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
                        MTF_ts(end) = Durations(2);
                        
                        
                        %% 2) GET DATA
                        
                        % Preallocate
                        stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            FR_trs = nan(numel(trans_TDidx),numel(MTF_ts));
                            sp01 = zeros(numel(trans_TDidx),Durations(2));
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                
                                psth(pstid==unique(TrialData(pst_TDidx,:).trID)',:,it) = ...
                                    Stream_FRsmooth( t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );
                                
                                stim(pstid==unique(TrialData(pst_TDidx,:).trID)',:,it) = ...
                                    SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );                                
                                
                                sp=[]; sp = spiketimes( spiketimes>=t1(trans_TDidx(it)) ...
                                    & spiketimes<t3(trans_TDidx(it)) ) - t2(trans_TDidx(it)) - 1;
                                
                                sp01(it,sp(sp>0)+1) = 1;
                                
                                for ii = 1:numel(MTF_ts)
                                    FR_trs(it,ii) = sum(sp>=0 & sp<=(MTF_ts(ii)),2)/(MTF_ts(ii)/1000);
                                end
                                
                            end %it
                                                        
                            if it>=minTrs
                                % Store FR data
                                for ii = 1:numel(MTF_ts)
                                    eval( sprintf(' %s_data(istim,pstid,1) = mean(FR_trs(:,ii),1);         ',MTF_ts_str{ii}) )
                                    eval( sprintf(' %s_data(istim,pstid,2) = median(FR_trs(:,ii),1);       ',MTF_ts_str{ii}) )
                                    eval( sprintf(' %s_data(istim,pstid,3) = std(FR_trs(:,ii),1)/sqrt(it); ',MTF_ts_str{ii}) )
                                end
                                % Store sp01 data for VS
                                BinarySpikeData{istim,pstid} = sp01;
                            end
                            
                            
                        end %pstid (prev stim id)
                        
                        
                    end %ITIflag (Trial or ITI)                    
                end %istim
                
                
                %% Save this unit's FR data into the data structure
                
                for ii = 1:numel(MTF_ts)
                    FR_Population(N).(MTF_ts_str{ii}) = eval([MTF_ts_str{ii} '_data']);
                end
                
                %% Caluclate VS for this unit 
                
                VSdata_pdc2pdc = get_VS_periodic(BinarySpikeData,AMrates);
                
                
                % For IR stimuli
                VSdata_pdc2ir = get_VS_irregular(BinarySpikeData,AMrates);
                
                
                
                
                %% Plot this unit's context-separated MTF
                
                if PLOT_EACH_MTF
                
                which_ts = 4;
                
                thisData = FR_Population(N).(MTF_ts_str{which_ts});
                
                %%            TO PERIODIC  
                
                hfp = figure;
                set(hfp,'Position',fullscreen)
                
                %   * * * * * * * *     PERIODIC     * * * * * * * * 
                % FR
                subplot(2,3,1)
                hold on
                for ipst = [3 6]
                    errorbar(2:6, thisData(2:6,ipst,2)',thisData(2:6,ipst,3)',...
                        '--','Color',colors(ipst,:),'LineWidth',2);
                    ip(ipst/3)=plot(2:6,thisData(2:6,ipst,2),...
                        'o','Color',colors(ipst,:),'LineWidth',2,'MarkerSize',15);
                    errorbar(1, thisData(1,ipst,2)',thisData(1,ipst,3)',...
                        '--','Color',colors(ipst,:),'LineWidth',2);
                    plot(1,thisData(1,ipst,2),...
                        'o','Color',colors(ipst,:),'LineWidth',2,'MarkerSize',15);
                end
                set(gca,'XTick',1:6,'XTickLabel',{'Unmod.' AMrates},'xlim',[0 7])
                ylim([0 1+ceil(max(max(thisData(:,:,2)+thisData(:,:,3))))])
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('FR response (sp/s)')
                title('Periodic rMTF')
                legend(ip,{'following  4 Hz' 'following 32 Hz'},'Location','best')
                hold off; clear ip
                
                % RESPONSE DIFFERENCES
                subplot(2,3,2)
                hold on
                plot([0 7],[0 0],'-k','LineWidth',1)
                errorbar(2:6, thisData(2:6,6,2)'-thisData(2:6,3,2)', mean(thisData(2:6,[3 6],3),2)',...
                    '--k','LineWidth',2)
                ip(1)=plot(2:6,thisData(2:6,6,2)-thisData(2:6,3,2),...
                    'ok','LineWidth',2,'MarkerSize',15);
                errorbar([3 6], thisData([3 6],8,2)'-thisData([3 6],7,2)',mean(thisData([3 6],[7 8],3),2)',...
                    '--','Color',0.6*ones(1,3),'LineWidth',2);
                ip(2)=plot([3 6],thisData([3 6],8,2) - thisData([3 6],7,2),...
                        'o','Color',0.6*ones(1,3),'LineWidth',2,'MarkerSize',15);
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[0 7])
                ylim([-10 10])
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('Difference of FR responses')
                title('diff rMTF (32hz-4hz)')
                legend(ip,{'pdc-pdc (32hz-4hz)' 'IR-pdc (DB-AC)'},'Location','best')
                
                
                % VS
                
                % MTF
                subplot(2,3,4)
                hold on
                for ipst = [3 6]
                    ip(ipst/3)=plot(2:6,[VSdata_pdc2pdc(1:5,ipst/3).VS],...
                        'o--','Color',colors(ipst,:),'LineWidth',2,'MarkerSize',15);
                end
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[0 7])
                ylim([0 1])
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('VS')
                title('Periodic vsMTF')
                legend(ip,{'following  4 Hz' 'following 32 Hz'},'Location','best')
                hold off; clear ip
                % Diff
                subplot(2,3,5)
                hold on
                plot([0 7],[0 0],'-k','LineWidth',1)
                plot(2:6,[VSdata_pdc2pdc(1:5,2).VS] - [VSdata_pdc2pdc(1:5,1).VS],...
                    'ok--','LineWidth',2,'MarkerSize',15);
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[0 7])
                ylim([-1 1].*0.3)
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('Difference of VS responses')
                title('diff vsMTF (32hz-4hz)')
                hold off; clear ip
                
                
                
                %   * * * * * * * *     IRREGULAR     * * * * * * * * 
                
                % IR to periodic transitions (to ITI blocks)
                subplot(2,3,3)
                hold on
                prevStim = [7 8]; %currently pp_parse_data isn't grabbing ITI blocks following Warns
                for ipst = prevStim
                    errorbar(1:2, thisData([3 6],ipst,2)',thisData([3 6],ipst,3)',...
                        '--','Color',colors(ipst,:),'LineWidth',2);
                    ip(ipst==prevStim)=plot(1:2,thisData([3 6],ipst,2),...
                        'o','Color',colors(ipst,:),'LineWidth',2,'MarkerSize',15);
                end
                set(gca,'XTick',[1 2],'XTickLabel',AMrates([2 5]),'xlim',[0 3])
                ylim([0 1+ceil(max(max(thisData(:,:,2)+thisData(:,:,3))))])
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('FR response (sp/s)')
                title('IR to periodic transitions')
                legend(ip,{'following AC' 'following DB'},'Location','best')
                hold off; clear ip
                
                
                suptitle(sprintf('%s %s %i %i (%s) R=%i\n%s',subject,session,channel,clu,unType,RobustUn,MTF_ts_str{which_ts}))
                
                
                % SAVE FIGURE
                
                savedir = fullfile(fn.processed,'ContextMTFs');
%                 if RobustUn==0
%                     savedir = [savedir '/notRobust'];
%                 end
                if ~exist(savedir,'dir')
                    mkdir(savedir);
                end
                savename = sprintf('ContextMTFs_%s_%s_ch%i_%i_%s_%s',subject,session,channel,clu,unType,MTF_ts_str{which_ts});
                print_eps_kp(hfp,fullfile(savedir,savename))
                hold off
                
                
                
                
                %%            TO IRREGULAR
                
                hfi = figure;
                set(hfi,'Position',fullscreen)
                hold on
                
                
                % Periodic to IR transitions
                PdcWeightedAvg = sum( mean(thisData(2:6,[3 6],2),2,'omitnan')' .* ((1./AMrates)/sum(1./AMrates)) );
                
                % FR
                % plotting as in Bartlett
                subplot(2,3,1)
                hold on
                ip(1)=plot([1 7],[PdcWeightedAvg PdcWeightedAvg],':k','LineWidth',0.5);
                for ist = 7:8
                    errorbar(2:6, thisData(ist,2:6,2)',thisData(ist,2:6,3)',...
                        '--','Color',colors(ist,:),'LineWidth',2)
                    ip(ist-5)=plot(2:6,thisData(ist,2:6,2),...
                        'o','Color',colors(ist,:),'LineWidth',2,'MarkerSize',15);
                end
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[1 7])
                ylim([0 1+ceil(max(max(thisData(:,:,2)+thisData(:,:,3))))])
                xlabel('AMrate (Hz) of PRECEDING STIMULUS')
                ylabel('FR response (sp/s), split by preceding stim')
                legend(ip,{'pdc linear prediciton' 'AC' 'DB'},'Location','best')
                title('IR Block responses')
                hold off; clear ip
                
                
                % VS
                VSdata_pdc2ir
                
                % MTF
                subplot(2,3,4)
                hold on
                for ipst = 2:6
                    ip(ipst-1)=plot(ipst,mean([VSdata_pdc2ir(1:2,ipst-1).VS],'omitnan'),...
                        'o--','Color',colors(ipst,:),'LineWidth',2,'MarkerSize',15);
                end
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[0 7])
                ylim([0 1])
                xlabel('AMrate (Hz) of PRECEDING STIMULUS')
                ylabel('VS')
                title('Irregular vsMTF (avg across pds)')
                legend(ip,{'following  2 Hz' 'following  4 Hz' 'following  8 Hz' 'following  16 Hz' 'following 32 Hz'},'Location','best')
                hold off; clear ip
                % Diff
                subplot(2,3,5)
                hold on
                plot([0 7],[0 0],'-k','LineWidth',1)
                plot(2:6,[VSdata_pdc2pdc(1:5,2).VS] - [VSdata_pdc2pdc(1:5,1).VS],...
                    'ok--','LineWidth',2,'MarkerSize',15);
                set(gca,'XTick',2:6,'XTickLabel',AMrates,'xlim',[0 7])
                ylim([-1 1].*0.3)
                xlabel('AMrate (Hz) of THIS STIMULUS')
                ylabel('Difference of VS responses')
                title('diff vsMTF (32hz-4hz)')
                hold off; clear ip
                
                
                
                
                end
                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects

if nargin<0
    return
end




end %function




