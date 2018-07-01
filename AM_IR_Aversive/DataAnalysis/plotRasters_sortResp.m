function plotRasters_sortResp(optargin)
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

close all

%!!!!!!!!!!!
SUonly = 1;
%!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;


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



%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];

for ist=1:7
    
    hf(ist) = figure;
    if nargin>0
        set(gcf,'Position',largerect)
    else
        set(gcf,'Position',halfscreen)
    end
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xtick',[],'ytick',[])
    hs(ist,2) = subplot(9,1,2:6);
    set(gca,'xtick',[],'ytick',[])
    hs(ist,3) = subplot(9,1,7:9);
    
end

add_y = zeros(1,7);
N_un  = zeros(1,7);

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


%%

% Remove non-responsive units
[sigUnits,Resp] = identifyResponsiveUnits(Resp);
Resp = Resp(sigUnits);

yyy = cellfun(@(xxx) mean(xxx,1,'omitnan'), {Resp.FR_raw_tr},'UniformOUtput',false);
Resp_meanFRs = vertcat(yyy{:});


% Plot one stimulus at a time, without separating previous stimuli

for stid = 2:8
    
    [~,isrt] = sort(Resp_meanFRs(:,stid));
    Resp_sort = Resp(isrt);
    
    GroupPSTH{stid-1} = nan(numel(Resp_sort),3000);
    
    if nargin>0
        theseUnits = find(strcmp({Resp_sort.Session},'QA') & [Resp_sort.Channel]==3);
    else
        theseUnits = 1:numel(Resp_sort);
    end
%     theseUnits = theseUnits(1);
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
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        
        % Get all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
        allStim = unique(TrialData.trID(all_TDidx));
        
        if ~ismember(stid,allStim), continue, end
        
        
        %% FIRST, COLLECT AND SET SOME STIMULUS INFO
        
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
                
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        t3 = t2 + Duration;
        
        % Preallocate
        raster_x = []; 
        raster_y = [];
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
            
        end %it
        
        
        % Skip plotting if too few trials
        if it<minTrs, continue, end
        
        
        %% Add to plots
        
        figure(hf(stid-1)); hold on
        
        % Stimulus
        subplot(hs(stid-1,1)); hold on
        plot(0:Duration, mean(stim,1,'omitnan'),...
            'LineWidth',1)
        hold off
        
        % Raster
        if numel(raster_x)>0
            
        subplot(hs(stid-1,2)); hold on
            plot(raster_x, raster_y + add_y(stid-1),...
                '.','MarkerSize',5) %,'Color',pstcolors(pst_idx,:)
            hold off
            
            add_y(stid-1) = add_y(stid-1) + max(raster_y) ;
        end
        
        % PSTH
        subplot(hs(stid-1,3)); hold on
        plot(0:Duration, mean(psth,1,'omitnan'),...
            'LineWidth',1);
        hold off
        
        N_un(stid-1) = N_un(stid-1) + 1;
        GroupPSTH{stid-1}(N_un(stid-1),1:Duration+1) = mean(psth,1,'omitnan');
        
        
        
    end %iUn
    
    
    % Trim grand psth to proper duration
    GroupPSTH{stid-1} = GroupPSTH{stid-1}(1:N_un(stid-1),1:Duration+1);
    
    % Finish plots
    subplot(hs(stid-1,3)); hold on
    plot(0:Duration,mean(GroupPSTH{stid-1},1),'k','LineWidth',3);
    ylim([0 80])
    xlabel('Time (ms)')
    ylabel('Spikes/s')
    hold off
    
    subplot(hs(stid-1,2)); hold on
    ylim([0 add_y(stid-1)+1]);
    ylabel([num2str(N_un(stid-1)) ' units'])
    hold off
    
    linkaxes(hs(stid-1,:),'x')
    xlim([0 Duration])
    
end % istimID



%% Save figure

savedir = fullfile(fn.processed,'Rasters_Grouped');
if ~exist(savedir,'dir')
    mkdir(savedir)
    mkdir([savedir '/eps'])
    mkdir([savedir '/svg'])
end

for ist = 1:7
    figure(hf(ist)); hold on
    suptitle(sprintf('%s, all responsive SU',Info.stim_ID_key{ist+1}));
    
    subplot(hs(ist,3)); hold on
    
    if nargin>0
        ylim([0 40])
        savename = sprintf('%s_%i_%i_%s',session,channel,clu,Info.stim_ID_key{ist+1});
    else
        ylim([0 80])
        savename = sprintf('AllSU_%s',Info.stim_ID_key{ist+1});
        
        % Save Group PSTH data
        save(fullfile(savedir,'AllSU_GroupPSTH'),'GroupPSTH','-v7.3')
    end
        
    print_eps_kp(hf(ist),fullfile([savedir '/eps'],savename))
    print_svg_kp(hf(ist),fullfile([savedir '/svg'],savename))
end




%% Inspect the dynamic range (max-min FR) across stimuli

keyboard 

figure; hold on

for ist = 1:7
    
    subplot(1,2,1); hold on
    plot(ist*ones(size(GroupPSTH{ist},1)),range(GroupPSTH{ist},2),'ok')
    hold off
    
    subplot(1,2,2); hold on
    plot(ist,range(mean(GroupPSTH{ist},1)),'ok')
    hold off

    dynRange(ist) = mean(range(GroupPSTH{ist},2));
end

Prediction = sum((1000./AMrates)/sum(1000./AMrates).*dynRange(1:5));


keyboard

end