function CorrScore_Resp(optargin)
%
%  pp_plot_rasters(subject, session, [channel, clu] )
%    Plots a raster and psth for each stimulus, separating trials by
%    the preceding stimulus. 
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2018-04
%



winLen = 30;


%% Load Resp data table 

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
clear q


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
narrow     = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];


figwidthscales = [1 1 1 1 1 1 1.937 1.937];

add_y = zeros(1,7);
N_un  = zeros(1,7);

% Set colors
colors = [   0 200 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        
alphval = 0.6;
dotsize = 12;


%%

minTrs   =  10;

% Sort Resp by average FR across stim
yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {Resp.FR_raw_tr},'UniformOutput',false);
Resp_meanFRs = mean(vertcat(yyy{:}),2,'omitnan');
[~,isrt] = sort(Resp_meanFRs);
Resp_meanFRs = Resp_meanFRs(isrt);
Resp_sort = Resp(isrt);

% Identify units labeled AM responsive
[sigUnits,Resp_sort] = identifyResponsiveUnits(Resp_sort);

if nargin>0
    theseUnits = find(strcmp({Resp_sort.Session},'QA') & [Resp_sort.Channel]==3 & [Resp_sort.Clu]==23);
    %     theseUnits = theseUnits(1);
else
    keyboard
    theseUnits = 1:17;%numel(Resp_sort);
end

for iUn = theseUnits
    
    % Get this unit's info
    subject = Resp_sort(iUn).Subject;
    session = Resp_sort(iUn).Session;
    channel = Resp_sort(iUn).Channel;
    clu     = Resp_sort(iUn).Clu;
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    if (iUn>1 && ~( strcmp(subject,Resp_sort(iUn-1).Subject) && strcmp(session,Resp_sort(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,Resp_sort(iUn-1).Subject) && strcmp(session,Resp_sort(iUn-1).Session) && channel==Resp_sort(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % Get stimulus params
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    % Get spiketimes
    spikes = Spikes.sorted(channel);
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    
    
    % Convert FR to z-score
    bs_smth = 20;
    Stream_FRsmooth = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
    
    
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    allStim = unique(TrialData.trID(all_TDidx))';
    
    STIMidx = find(allStim==2 | allStim==7 | allStim==8);
    
    
    close all
    % Set up figure
    hf = figure;
    set(gcf,'Position',largerect)
    hold on
    hs(1) = subplot(2,1,1);
    set(gca,'xlim',[0 9])
    hs(2) = subplot(2,1,2); hold on
    set(gca,'xlim',[0 9],'ylim',[-0.3 0.7])
    plot([0 9],[0 0],'k-'); hold off
    
    
    ymaxval = 0;
    score = nan(numel(allStim),2);
    
    for ist = 1:numel(allStim)
        
        stid = allStim(ist);
        
        
        %% Collect trial indices and timestamps
        
        if stid==3 || stid==6
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
            % Find Pdc trials that follow same rate during ITI
            TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
            
            TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)==1);
            TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
            
        else
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
            TDidx_iti = [];
        end
        
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        
        % Add ITI trials (shortened to match duration)
        if ~isempty(TDidx_iti)
            t2 = [t2; TrialData.onset(TDidx_iti)];
        end
        
        t3 = t2 + Duration;
        
        TDidx = [TDidx; TDidx_iti];
        
        
        %%
        % Preallocate
        raster_x = []; 
        raster_y = [];
        raster01 = zeros( numel(TDidx), Duration ); 
        stim   = nan( numel(TDidx), Duration+1 ); 
        psth   = nan( numel(TDidx), Duration+1 ); 
        
        
        % Collect spikes/FR/rms for this stimulus/unit
        for it = 1:numel(TDidx)
            
            psth(it,:) = ...
                Stream_FRsmooth( t2(it) : t3(it) );
            
            stim(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
            sp=[]; sp = spiketimes( spiketimes>t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) - 1;
            
            raster_x = [raster_x sp];
            raster_y = [raster_y it*ones(1,numel(sp))];
            
            if ~isempty(sp)
                raster01(it,sp+1) = 1;
            end
            
        end %it
        
        
        % Skip plotting if too few trials
        if it<minTrs, continue, end
        
        
        
        %% Calculate CorrScore
        
%         score(ist,1) = calculateCorrScore( raster01(:,1:floor(size(raster01,2)/2)), winLen, Info.stim_ID_key{stid}  ,0 );
%         score(ist,2) = calculateCorrScore( raster01(:,ceil(size(raster01,2)/2):end), winLen, Info.stim_ID_key{stid} ,0 );
        
        score(ist,1) = calculateCorrScore( raster01, winLen, Info.stim_ID_key{stid} ,1 );
        
        
        %% Add to plots
        
        figure(hf); hold on
        
        % AM rate tuning curve
        subplot(hs(1)); hold on
        plot(stid, mean(mean(psth,1,'omitnan')),...
            'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',4)
%         plot(stid, mean(mean(psth(:,1:floor(size(raster01,2)/2)),1,'omitnan')),...
%             'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',4)
%         plot(stid, mean(mean(psth(:,ceil(size(raster01,2)/2):end),1,'omitnan')),...
%             'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',2)
        ylabel('Spikes/sec')
        hold off
        
        % Corr Score
        subplot(hs(2));  hold on
        plot(stid, score(ist,1),...
            'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',4)
%         plot(stid, score(ist,1),...
%             'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',4)
%         plot(stid, score(ist,2),...
%             'o','MarkerSize',dotsize,'MarkerFaceColor','none','MarkerEdgeColor',colors(stid,:),'LineWidth',2)
        ylabel('CorrScore')
        hold off
        
        
        suptitle(sprintf('%s   |   %s %s ch%i clu%i\nN=%i, avg FR = %0.2f Hz  |  Resp = %i  |  BMF = %i Hz ',...
            Info.stim_ID_key{stid}, subject,session,channel,clu,iUn,Resp_meanFRs(iUn),ismember(iUn,sigUnits),AMrates(Resp_sort(iUn).iBMF_FR) ))
        
        
        ymaxval = max(ymaxval,mean(mean(psth,1,'omitnan')));
        
        
        
    end %ist
    
    
    %% Save figure
    
    subplot(hs(1)); hold on
    ylim([0 ymaxval*1.25])
    
    subplot(hs(2)); hold on
    ylim([min(min(score))-0.05 max(max(score))+0.05])
    
    savedir = fullfile(fn.processed,'CorrScore_tmp',['win' num2str(winLen)]);
    if ~exist(savedir,'dir')
        mkdir(savedir)
%         mkdir([savedir '/eps'])
%         mkdir([savedir '/svg'])
    end
    
    savename = sprintf('un%i_%s_%s_%i_%i',iUn,subject,session,channel,clu);
    
    print_eps_kp(hf,fullfile(savedir,savename))
    
%     print_eps_kp(gcf,fullfile(savedir,[savename '_DB']))
    
    
end % iUn


keyboard 






end