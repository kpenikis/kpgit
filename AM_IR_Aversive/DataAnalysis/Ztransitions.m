function Ztransitions(SUBJECT, SESSION, select_clu)
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

global t_win wincolors pstcolors figsize_psth 

%!!!!!!!!!!!!!!!!!
PLOT_PSTH = 0;
%!!!!!!!!!!!!!!!!!
nIterations = 1000;
PLOT_RND_PSTH = 0;
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
FRcutoff = 1;
%!!!!!!!!!!!!!!!!!
tstep = 100;
t_win = 1:tstep:1001;


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Filter Unit files to just this session and sort by baseline FR
UnitData = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
[~, baseFRrank] = sort([UnitInfo.BaseFR]);
UnitInfo = UnitInfo(baseFRrank,:);
UnitData = UnitData(baseFRrank);

% Label as RS or NS (SpkShape)
UnitInfo = labelRSNS(UnitInfo);


% Also find trials to skip 
Channels = unique([UnitInfo.Channel]);
artifactTrs = [];
for ich = Channels'
    artifactTrs = [artifactTrs Info.artifact(ich).trials'];
end
artifactTrs = unique(artifactTrs);


% Get raster data
[UnitInfo, UnitData, Info, TrialData, Clusters, StimResp ] = collectRasterDataSession(SUBJECT,SESSION);

rng('shuffle')


%% 
% Set up figure options

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
figsize_psth = [1 scrsz(4) 2*scrsz(3)/3 scrsz(4)];
figsize1 = [1 1 scrsz(3) scrsz(4)/3];

colors = [  0   0   0;...
           84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [colors; 0.7.*bone(2)];

wincolors = flipud(winter(numel(t_win)-1));

pstcolors = [0 0 0; 0.7 0.7 0.7];



%% 
% Set up results table

Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.unN      = 0;
Zdata.spl      = 0;
Zdata.lpn      = 0;
Zdata.stim     = 0;

N = 0;



for iUn = 1:numel(UnitData)
    
    %%% still must edit for merged units
    %     if numel(UnitInfo(iUn,:).Session{:})>2  %strncmp(UnitInfo.RespType{iUn},'merged',6)
    %         continue
    %     end
    
    subject = UnitInfo(iUn,:).Subject{:};
    session = UnitInfo(iUn,:).Session{:};
    channel = UnitData(iUn).Channel(1);
    clu     = UnitData(iUn).Clu(1);
    
    % Load data files
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % Get spiketimes and shift based on calculated integration time (KS)
    iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
    Shank = Clusters(iClu).shank;
    spiketimes = round(Clusters(iClu).spikeTimes*1000)';
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~
    % Convert FR to z-score
    
    % Check the duration of silence at the beginning
    if ((TrialData.offset(1) - TrialData.onset(1))/1000) < 15
        keyboard
    end
    
    bs_smth = 20;
    [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    
    % Check if unit has very low overall firing rate
    if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
        disp(' few spiking events...skipping\n')
        continue
    else
        fprintf('~~~~ ch %i clu %i\n',channel,clu)
    end
    
    N = N+1;
    
    
    %%
    
    % Get sound parameters
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    [all_TDidx,Ntrials1,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % For each STIMULUS
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    for stid = unique(TrialData.trID(all_TDidx))'
        
        
        %% First, get stimulus indices and timestamps
        
        st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==stid);
        
        %%%  plot for Trial and ITI stimuli separately
        ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
        
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
            clear t1 t2 t3 Durations
            t2 = TrialData.onset(st_TDidx);
            t3 = TrialData.offset(st_TDidx);
            Durations(2) = mode(diff([t2 t3],1,2));
            t3 = t2 + Durations(2);
            Durations(1) = mode(diff([(TrialData.onset(pst_TDidx)) t2],1,2)); %doesn't make much sense when prev stim is ITI
            t1 = t2 - Durations(1);
            
            
            %% GET DATA 
            
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
                
                % Make raster vectors for this transition
                eval( sprintf('raster_x_%i = [];',pst_idx) )
                eval( sprintf('raster_y_%i = [];',pst_idx) )
                
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
                    
                    eval( sprintf('raster_x_%i = [raster_x_%i sp];',pst_idx,pst_idx) )
                    eval( sprintf('raster_y_%i = [raster_y_%i it*ones(1,numel(sp))];',pst_idx,pst_idx) )
                    
                end %it
                
                % Save N trials
                Ntrials(1,pst_idx) = it;
                
                % Save the psth data for each time window
                for iw = 2:numel(t_win)
                    IR_psth_wins(pst_idx,:,iw-1) = mean( psth( pst_idx, (t_win(iw-1):t_win(iw)-1)+Durations(1), :) ,3,'omitnan');
                end
                
                
            end %pstid (prev stim id)
            
            
            % Save the difference in spikes between conditions
            TestDiffSpikeCount = nan(1,numel(t_win)-1);
            if size(IR_psth_wins,1)>1
                for iw = 2:numel(t_win)
                    % Get cumulative count of difference in n
                    % spikes
                    TestDiffSpikeCount(1,iw-1) = sum(abs( diff(IR_psth_wins([1 size(IR_psth_wins,1)],:,iw-1)./1000,1) )); %if more than 2 prev stim, just save the first and last prev stim conditions
                end
            end
            
            
            
            %% Plot raster/PSTH around transition
            
            if PLOT_PSTH
                
                stimstring  = sprintf('Transition to %s stimulus',Info.stim_ID_key{stid});
                TitleString = sprintf('%s\n%s %s  |  shank %i ch %i clu %i (%s)',...
                    stimstring, subject,session,Shank,channel,clu,UnitInfo.SpkShape{iUn});
                
                for pstid = unique(TrialData(pst_TDidx,:).trID)'
                    ContextString{pstid==unique(TrialData(pst_TDidx,:).trID)'} = Info.stim_ID_key{pstid};
                end
                if size(IR_psth_wins,1)>2
                    ContextString = ContextString([2 end]);
                    % NOTE HERE, overwriting data!
                    IR_psth_wins  = IR_psth_wins([2 end],:,:);
                    psth          = psth([2 end],:,:);
                    stim          = stim([2 end],:,:);
                    raster_x_1 = raster_x_2; 
                    raster_y_1 = raster_y_2; 
                    raster_x_2 = raster_x_5; 
                    raster_y_2 = raster_y_5;
                end
                
                plotTransitionPSTH( IR_psth_wins, psth, stim, ContextString, TitleString, Durations, ymaxval, raster_x_1, raster_x_2, raster_y_1, raster_y_2 )
                
            end %PLOT_PSTH
            
            
            
            
            %%
            %++++++++++++++++++++++++++++++++++++++++++++++++++
            %  Bootstrap random trials, for comparison
            %++++++++++++++++++++++++++++++++++++++++++++++++++
            
            % Skip stimuli that don't have at least 10 trials,
            % or that don't have exactly 2 preceding stimulus
            % types, for now
            if (numel(unique(TrialData(pst_TDidx,:).trID))==2) && all(Ntrials>=minTrs)
                                
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
%                             IR_psth_wins(ihalf,:,iw-1) = mean( psth(ihalf,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
                            IR_psth_wins(ihalf,:,iw-1) = mean( psth(ihalf,(t_win(iw-1):t_win(iw)-1)+Durations(1), :) ,3,'omitnan');
                        end %iw
                        
                        
                    end %ihalf
                    
                    
                    % Calculate the difference between the PSTHs for
                    % each time window
                    for iw = 2:numel(t_win)
                        DistributionDiffSpikeCount(iteration,iw-1) = sum(abs( diff(IR_psth_wins(:,:,iw-1)./1000,1) ));
                    end
                    
                end %iteration
                
                
                %%
                %++++++++++++
                %  Z SCORE
                %++++++++++++
                
                % Calculate how far the Test spike counts are from
                % the bootstrapped distributions
                zScores = (TestDiffSpikeCount - mean(DistributionDiffSpikeCount,1)) ./ std(DistributionDiffSpikeCount,1);
                
                
                % Save result to table
                
                if ~any(strncmpi(Zdata.Properties.VariableNames,'zscore',6))
                    for iw = 2:numel(t_win)
                        Zdata.(sprintf('zScore_%i',t_win(iw)-1)) = 0;
                    end
                end
                
                Zdata_addrow = {subject session channel clu N dBSPL LP stid };
                for iw = 2:numel(t_win)
                    Zdata_addrow{end+1} = zScores(iw-1);
                end
                
                Zdata = [Zdata; Zdata_addrow];
                
                
                clear zScores
                
                
            end %stimulus filter for zscore analysis
            
            
            
            %% Plot last iteration of random split, if directed
            
            if PLOT_RND_PSTH && (numel(unique(TrialData(pst_TDidx,:).trID))==2) && all(Ntrials>=minTrs)
                
                stimstring  = sprintf('Transition to %s stimulus - RANDOM SPLIT',Info.stim_ID_key{stid});
                TitleString = sprintf('%s\n%s %s  |  shank %i ch %i clu %i (%s)',...
                    stimstring, subject,session,Shank,channel,clu,UnitInfo.SpkShape{iUn});
                ContextString = {'Random half 1' 'Random half 2'};
                
                plotTransitionPSTH( IR_psth_wins, psth, stim, ContextString, TitleString, Durations, ymaxval, raster_x_1, raster_x_2, raster_y_1, raster_y_2 )
                
            end %PLOT_RND
            
            
            
        end %ITIflag (Trial or ITI)
    end %stid (this stim id)
    
end %iUn


keyboard

% Save Zdata table
Zdata(1,:) = [];

% if SUonly==1
%     savedir = fullfile(fn.processed,'Transitions','SU');
% else
%     savedir = fullfile(fn.processed,'Transitions','allU');
% end
% if ~exist(savedir,'dir')
%     mkdir(savedir)
% end
% 
% writetable(Zdata,fullfile(savedir,sprintf('PdcTransitions_zSpkCtDiffs_%imswin_RobustUn',tstep)));


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
savename = sprintf('zSpkCtDiffs_%s_%imswin_individualUnits',SESSION,tstep);
print_eps_kp(hf,fullfile(fn.processed,'Transitions',savename))


keyboard

[~,idx] = sort(Zdata.zScore_100);
Zdata(idx,:)


end %function




