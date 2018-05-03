function plotAllRasters( )
%
%  plotAllRasters( ) 
%    Plots a raster and psth for each stimulus, separating trials by
%    the preceding stimulus. 
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2018-05
%

close all

%!!!!!!!!!!!
SUonly = 1;
%!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;


%% Load Resp data table (with RCorr results)

fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'RespStruct_allSU'));
Resp = q.Resp;

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;
AMrates = [2 4 8 16 32];
IRstr = {'AC' 'DB'};

N=0;


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
narrow     = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];

for iir = 1:2
for ir=1:5
    
    hf(iir,ir) = figure;
    set(gcf,'Position',fullscreen,'colormap',prism(40))
    hold on
    hs(iir,ir,1) = subplot(9,1,1);
    set(gca,'xtick',[],'ytick',[])
    hs(iir,ir,2) = subplot(9,1,2:6);
    set(gca,'xtick',[],'ytick',[])
    hs(iir,ir,3) = subplot(9,1,7:9);
    
    suptitle([num2str(AMrates(ir)) ' Hz - IR ' IRstr{iir}])
    add_y(iir,ir) = 0;
    
end
end


% Set colors
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        
% pstcolors = colors(unique(TrialData(pst_TDidx,:).trID),:);
% wincolors = flipud(winter(numel(t_win)-1));

alphval = 0.6;
dotsize = 120;




%%  SUBJECTS

subjects = {'WWWlf_253395' 'WWWf_253400' };

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

SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));

Sessions = [];
for ifn = 1:numel(SpkFns)
    if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
        Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
    end
end
% Sessions = flipud(Sessions);


% Step through each session
for sess = Sessions'
    
session = char(sess);

SessResp = Resp(strcmp({Resp.Subject},subject)&strcmp({Resp.Session},session));
if isempty(SessResp),  continue,   end

% Here add call to a function that outputs the indices of SessResp
% that contain a responsive unit. With this method, the thresholds can
% be flexibly changed.
[sig_ids,SessResp] = identifyResponsiveUnits(SessResp);

Channels = unique([SessResp(sig_ids).Channel]);

if isempty(Channels), continue, end



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


%% STEP THROUGH EACH CHANNEL

for channel = Channels
        
    % Find clus to plot
    Clus = SessResp([SessResp.Channel]==channel).Clu;
    
    spikes = Spikes.sorted(channel);
    
    
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
        
        ThisResp = SessResp([SessResp.Channel]==channel & [SessResp.Clu]==clu);
        
        
        
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
                
                all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                
                allStim = unique(TrialData.trID(all_TDidx));
                IRstim = allStim(allStim>6);
                
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                for stid = IRstim'
                    
                    iir = stid-6;
                    
                    %% FIRST, COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==stid);
                    
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
                        if Durations(1)>1005
                            keyboard
                        end
                        t1 = t2 - Durations(1);
                        tstep = 100;
                        t_win = 1:tstep:1001;
                        t_FF = [-1000:10:999]';
                        
                        
                        % Preallocate
%                         legstr = cell(1,numel(unique(TrialData(pst_TDidx,:).trID))); clear ip
                        vars = whos;
                        cellfun(@clear,({vars(~cellfun(@isempty,regexp({vars.name},'raster_*'))).name}))
                        stim = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        psth = nan( numel(unique(TrialData(pst_TDidx,:).trID)), sum(Durations)+1, 100 ); % ( previous stim, duration, trials )
                        FF_mean  = nan(numel(unique(TrialData(pst_TDidx,:).trID)),numel(t_FF));
                        FF_var   = nan(numel(unique(TrialData(pst_TDidx,:).trID)),numel(t_FF));
                        
                        
                        %% GET DATA AND PLOT
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            pst_idx = find(pstid==unique(TrialData(pst_TDidx,:).trID)');
                            
                            % Make raster vectors for this transition
                            eval( sprintf('raster_x_%i = [];',pst_idx) )
                            eval( sprintf('raster_y_%i = [];',pst_idx) )
                            
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            FF_bincounts = nan(numel(trans_TDidx),numel(t_FF));
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                                                
                                psth(pst_idx,:,it) = ...
                                    Stream_FRsmooth( t1(trans_TDidx(it)) : t3(trans_TDidx(it)) );
                                
                                stim(pst_idx,:,it) = ...
                                    SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) )...
                                    ./ max(SoundStream(1, t1(trans_TDidx(it)) : t3(trans_TDidx(it)) ));                                
                                
                                sp=[]; sp = spiketimes( spiketimes>=t1(trans_TDidx(it)) ...
                                    & spiketimes<t3(trans_TDidx(it)) ) - t2(trans_TDidx(it)) - 1;
                                
                                eval( sprintf('raster_x_%i = [raster_x_%i sp];',pst_idx,pst_idx) )
                                eval( sprintf('raster_y_%i = [raster_y_%i it*ones(1,numel(sp))];',pst_idx,pst_idx) )
                                
                                FF_bincounts(it,:) = sum(sp>=t_FF & sp<(t_FF+tstep),2)';
                                
                                
                            end %it
                            
                            % Collect FF data
                            FF_mean(pst_idx,:) = mean(FF_bincounts,1);
                            FF_var(pst_idx,:)  = var(FF_bincounts,1);
                            
                            
                            % Skip plotting if too few trials
                            if it<minTrs
                                continue
                            end
                            
                            
                            %% Add to plots
                            
                            figure(hf(iir,pstid-1)); hold on                            
                            
                                                        
                            % Stimulus
                            subplot(hs(iir,pstid-1,1)); hold on
                            plot(-Durations(1):Durations(2), mean(stim(pst_idx,:,:),3,'omitnan'),...
                                'LineWidth',1)
                            hold off
                            
                            % Raster
                            if numel(eval(['raster_y_' num2str(pst_idx)]) )>0
                                
                                subplot(hs(iir,pstid-1,2)); hold on
                                plot(eval(['raster_x_' num2str(pst_idx)]),eval(['raster_y_' num2str(pst_idx)]) + add_y(iir,pstid-1),...
                                    '.','MarkerSize',5) %,'Color',pstcolors(pst_idx,:)
                                hold off
                                
                                add_y(iir,pstid-1) = add_y(iir,pstid-1)+ max( eval(['raster_y_' num2str(pst_idx)]) ) ;
                                
                            end
                            
                            % PSTH
                            subplot(hs(iir,pstid-1,3)); hold on
%                             ip(pst_idx) = 
                            plot( -Durations(1):Durations(2), mean(psth(pst_idx,:,:),3,'omitnan') ,...
                                'LineWidth',1);
                            hold off
%                             legstr{pst_idx} = [Info.stim_ID_key{pstid} ', n=' num2str(sum((TrialData(pst_TDidx,:).trID==pstid)))];
                                                        
                            
                            
                        end %pstid (prev stim id)
                    end %ITIflag (Trial or ITI)
                    
                end %stid (this stim id)

            end %lpn
        end %spl
        
    end %clu
end %channel
end %sess
end %subject


%% Finish plot settings


linkaxes(hs,'x')
xlim([-Durations(1) Durations(2)])

for iir = 1:2
    for ir=1:5
        figure(hf(iir,ir)); hold on
        
        subplot(hs(iir,ir,2)); hold on
        ylim([0 add_y(iir,ir)+1]);
        hold off
        
        subplot(hs(iir,ir,3)); hold on
        ylim([0 100])
        xlabel('Time from transition (ms)')
        ylabel('Spikes/s')
        
        hold off
        
    end
end


%% Save figure
keyboard
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



end