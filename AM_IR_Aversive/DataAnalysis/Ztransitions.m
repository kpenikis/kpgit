function Ztransitions(select_subject, select_session, select_channel, select_clu)
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
%!!!!!!!!!!!!!!!!!
RobustUn =  1;
%!!!!!!!!!!!
nIterations = 1000;
PLOT_RND = 0;
%!!!!!!!!!!!
FRcutoff =  2;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
tstep = 100;
t_win = 1:tstep:1001;

%% Load Resp data table (with RCorr results)

try
    fn = set_paths_directories('','',1);
    Resp = readtable(fullfile(fn.processed,'Resp_allSU'));
catch
    RobustUn = 0;
end


rng('shuffle')

%% 
% Set up figure options

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

wincolors = flipud(winter(numel(t_win)-1));


%% 
% Set up results table

Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.unType   = ' ';
Zdata.unN      = 0;
Zdata.spl      = 0;
Zdata.lpn      = 0;
Zdata.stim     = 0;

N = 0;



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

if nargin>2 && exist('select_channel','var')
    Channels = select_channel;
else % group data
    if RobustUn %skip sessions with no Robust SUs
        SessResp = Resp(strcmp(Resp.Subject,subject)&strcmp(Resp.Session,session),:);
        SessResp = SessResp(SessResp.AudResp_rc==1,:);
        if isempty(SessResp),  continue,   end
        Channels = unique(SessResp.ch)';
    else
        Channels =  [1:7 9:16];
    end
end


%%
% Load data files
fn = set_paths_directories(subject,session,1);
fprintf('Loading data for Session %s...\n',session)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(fn.processed,subject,filename));


% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);



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
        
        if isempty(spiketimes)
            continue
        end
        

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
        if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
%             disp(' few spiking events')
            continue
        else
            fprintf('ch %i clu %i\n',channel,clu)
        end
        
        N = N+1;
        
        
        %%
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
                [all_TDidx,Ntrials1,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                allStim = unique(TrialData.trID(all_TDidx));
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                for stid = 1:numel(allStim)
                    
                    
                    %% FIRST, COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==allStim(stid));
                    
                    %%%  plot for Trial and ITI stimuli separately
                    ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for is = 1%:numel(ITIflag)
                        
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
                        clear t1 t2 t3 Durations  
                        t2 = TrialData.onset(st_TDidx);
                        t3 = TrialData.offset(st_TDidx);
                        Durations(2) = mode(diff([t2 t3],1,2));
                        t3 = t2 + Durations(2);
                        Durations(1) = mode(diff([(TrialData.onset(pst_TDidx)) t2],1,2));
                        t1 = t2 - Durations(1);
                        
                        
                        %% GET DATA AND PLOT
                        
                        % Preallocate
                        legstr = cell(1,numel(unique(TrialData(pst_TDidx,:).trID))); clear ip
                        vars = whos;
                        cellfun(@clear,({vars(~cellfun(@isempty,regexp({vars.name},'raster_*'))).name}))
                        stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        IR_psth_wins = [];
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        Ntrials = nan(1,numel(unique(TrialData(pst_TDidx,:).trID)));
                        
                        jt=0;
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            pst_idx = find(pstid==unique(TrialData(pst_TDidx,:).trID)');
                            
%                             % Make raster vectors for this transition
%                             eval( sprintf('raster_x_%i = [];',pst_idx) )
%                             eval( sprintf('raster_y_%i = [];',pst_idx) )
                            
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                
                                jt=jt+1;
                                
                                psth(pst_idx,:,it) = ...
                                    Stream_FRsmooth( t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );
                                
                                stim(pst_idx,:,it) = ...
                                    SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );                                
                                
                                sp=[]; sp = spiketimes( spiketimes>=t1(trans_TDidx(it)) ...
                                    & spiketimes<t3(trans_TDidx(it)) ) - t2(trans_TDidx(it)) - 1;
                                
%                                 eval( sprintf('raster_x_%i = [raster_x_%i sp];',pst_idx,pst_idx) )
%                                 eval( sprintf('raster_y_%i = [raster_y_%i it*ones(1,numel(sp))];',pst_idx,pst_idx) )

                            end %it
                            
                            % Save N trials
                            Ntrials(1,pst_idx) = it;
                            
                            % Save the psth data for each time window
                            for iw = 2:numel(t_win)
                                IR_psth_wins(pst_idx,:,iw-1) = mean( psth( pst_idx, (t_win(iw-1):t_win(iw)-1)+Durations(1), :) ,3,'omitnan');
                            end
                            
                            
                        end %pstid (prev stim id)
                        
                        
                        % For now, plot the difference between the first 
                        % and last prev stim conditions 
                        TestDiffSpikeCount = nan(1,numel(t_win)-1);
                        if size(IR_psth_wins,1)>1
                            for iw = 2:numel(t_win)
                                % Get cumulative count of difference in n
                                % spikes
                                TestDiffSpikeCount(1,iw-1) = sum(abs( diff(IR_psth_wins([1 size(IR_psth_wins,1)],:,iw-1)./1000,1) ));
                            end
                        end
                        
                        
                        %++++++++++++++++++++++++++++++++++++++++++++++++++
                        %% FINISH Z-SCORE ANALYSIS
                        %++++++++++++++++++++++++++++++++++++++++++++++++++
                        
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
%                                 eval( sprintf('raster_x_%i = [];',ihalf) )
%                                 eval( sprintf('raster_y_%i = [];',ihalf) )
                                
                                thesetrs = randomtrs( 1:Ntrials(1,ihalf) );
                                randomtrs( 1:Ntrials(1,ihalf) ) = [];
                                
                                % Collect FR
                                for jj = 1:numel(thesetrs)
                                    
                                    psth(ihalf,:,jj) = Stream_FRsmooth(t1(thesetrs(jj)==pst_TDidx) : t3(thesetrs(jj)==pst_TDidx) );
                                    
                                    stim(ihalf,:,jj) = SoundStream( 1, t1(thesetrs(jj)==pst_TDidx) : t3(thesetrs(jj)==pst_TDidx) );
                                    
%                                     sp = []; 
%                                     sp = spiketimes( spiketimes>=t1(thesetrs(jj)==pst_TDidx) & spiketimes<t3(thesetrs(jj)==pst_TDidx) )...
                                        - t2(thesetrs(jj)==pst_TDidx) - 1;
%                                     eval( sprintf('raster_x_%i = [raster_x_%i sp];',ihalf,ihalf) )
%                                     eval( sprintf('raster_y_%i = [raster_y_%i jj*ones(1,numel(sp))];',ihalf,ihalf) )
%                                     
                                end
                                
                                for iw = 2:numel(t_win)
                                    % Save the psth data for each time window
%                                     IR_psth_wins(ihalf,:,iw-1) = mean(psth(ihalf,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
                                    IR_psth_wins(ihalf,:,iw-1) = mean( psth(ihalf,(t_win(iw-1):t_win(iw)-1)+Durations(1), :) ,3,'omitnan');
                                end %iw
                                
                                
                            end %ihalf
                            
                            
                            % Calculate the difference between the PSTHs for
                            % each time window
                            for iw = 2:numel(t_win)
                                DistributionDiffSpikeCount(iteration,iw-1) = sum(abs( diff(IR_psth_wins(:,:,iw-1)./1000,1) ));
                            end
                            
                        end %iteration
                        
                        
                        %% Z SCORE
                        
                        % Calculate how far the Test spike counts are from
                        % the bootstrapped distributions
                        zScores = (TestDiffSpikeCount - mean(DistributionDiffSpikeCount,1)) ./ std(DistributionDiffSpikeCount,1);
                        
                        
                        %% Save result to table
                        
                        if ~any(strncmpi(Zdata.Properties.VariableNames,'zscore',6))
                            for iw = 2:numel(t_win)
                                Zdata.(sprintf('zScore_%i',t_win(iw)-1)) = 0;
                            end
                        end
                        
                        Zdata_addrow = {subject session channel clu unType N spl lpn stid };
                        for iw = 2:numel(t_win)
                            Zdata_addrow{end+1} = zScores(iw-1);
                        end
                        
                        Zdata = [Zdata; Zdata_addrow];
                        
                        
                        clear zScores
                        
                        
                        end %stimulus filter for zscore analysis
                        
                        
                        
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
            end %lpn
        end %spl
    end %clu
end %channel
end %session
end %subject


% Save Zdata table
Zdata(1,:) = [];

if SUonly==1
    savedir = fullfile(fn.processed,'Transitions','SU');
else
    savedir = fullfile(fn.processed,'Transitions','allU');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

writetable(Zdata,fullfile(savedir,sprintf('PdcTransitions_zSpkCtDiffs_%imswin_RobustUn',tstep)));


% Plot the result
hf=figure;
plot([2.5 2.5],[0 1],'Color',[0.5 0.5 0.5])
hold on
plot([-2.5 -2.5],[0 1],'Color',[0.5 0.5 0.5])
for iw = numel(t_win):-1:2
    [p2,x2]=ecdf(Zdata.(sprintf('zScore_%i',t_win(iw)-1)));
    ip(iw-1)=plot(x2,p2,'LineWidth',5,'Color',wincolors(iw-1,:));
    legtext{iw-1} = sprintf('%i-%i ms',t_win(iw)-1-mode(diff(t_win)),t_win(iw)-1);
end

% title(subject)
xlabel('z-score of spike count difference when trials split by prev pdc rate')
ylabel('Probability')
legend(ip,legtext,'Location','southeast','Interpreter','none')
text(1,0.05,sprintf('N = %i',size(Zdata,1)))

% Save population responses figure
savename = sprintf('Population_PdcTransitions_zSpkCtDiffs_%imswin_RobustUn',tstep);
print_eps_kp(hf,fullfile(savedir,savename))


keyboard

[~,idx] = sort(Zdata.zScore_100);
Zdata(idx,:)


end %function




