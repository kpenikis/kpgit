function inspectRespTiming(SUBJECT,SESSION)
% 
% After meeting with Mike Long, looking at the data in various ways. Here,
% plot a few hand-picked sparse responding cells' responses to the 4 Hz 
% stimulus, sorted by resp phase. Plot responses to 2 and 8 Hz stimuli, as
% well as 4 Hz periods from aperiodic stimuli, in the same sequence. 
% 
% KP, 2019-03-07
% 


% close all
global fn AMrates rateVec_AC rateVec_DB 


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
UnitData = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
% [~, UIrank] = sort([UnitInfo.BaseFR]);
% UnitInfo = UnitInfo(UIrank,:);
% UnitData = UnitData(UIrank);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

% Sort units by FR (to help pick out units from allRastersSession plot) 
[~, inspks] = sort(cellfun(@length, {Clusters.spikeTimes}));
Clusters = Clusters(inspks);

UnitInfo = UnitInfo(inspks,:);
UnitData = UnitData(inspks);


% Sequence of response peaks, by eye, from Session raster
% %
cluIDs = [1650 1139 771 1230]; %68
% %


%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
smallsq  = [1 scrsz(4) scrsz(3)/2 scrsz(4)/2];

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
colors = colors.*0;


rasterdotsize  = 15;
% UnitColors = [0 0 0; jet(length(cluIDs)-1)];
UnitColors = [ 0 0 0; 0 0 1; 0 0.8 0.1; 0.9 0.1 0];

        
% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];

minTrs = 12;


%% For 4, 2, and 8 Hz stimuli

for stid = [3 2 4 5 6]
    
    % Prepare the figure
    hf(stid-1) = figure;
    set(hf(stid-1),'Position',smallsq,'NextPlot','add')
    hold on
    hs(1)=subplot(6,1,1);
    set(gca,'xtick',[],'ytick',[],'Color','none');
    box off
    hs(2)=subplot(6,1,2:6);
    set(gca,'Color','none','tickdir','out')
    box off
    
    
    % Preallocate
    N_un = 0;
    add_y = 0;
    
    
    for iClu = cluIDs
        
        iUn = [UnitData.Clu]==iClu;
        %%% still must edit for merged units
        %     if numel(UnitInfo(iUn,:).Session{:})>2  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        %         continue
        %     end
        
        channel = UnitData(iUn).Channel(1);
        clu     = UnitData(iUn).Clu(1);
        
        % Get spiketimes (KS)
        spiketimes = round(Clusters(([Clusters.maxChannel]==channel & [Clusters.clusterID]==clu)).spikeTimes*1000)';
        
        % Get sound parameters
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        fprintf(' analyzing ch %i clu %i\n',channel,clu)
        
        
        [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,20,'silence');
        
        %%
        % Find all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
        
        
        
        allStim = unique(TrialData.trID(all_TDidx));
        
        if sum(Ntrials < minTrs)==1
            keyboard
            all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
            allStim(Ntrials<minTrs)  = [];
            Ntrials(Ntrials<minTrs) = [];
        elseif  sum(Ntrials < minTrs)>1
            keyboard
        end
        
        %%
        % Get trial numbers for this stimulus
        
        if stid==3 || stid==6
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
            % Find Pdc trials that follow same rate during ITI
            %         TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
            %
            %         TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
            %         TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        else
            TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
            TDidx_iti = [];
        end
        
        
        % Timestamps of onsets and offsets
        
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        %
        %     % Add ITI trials (shortened to match duration)
        %     if ~isempty(TDidx_iti)
        %         t2 = [t2; TrialData.onset(TDidx_iti)];
        %         TDidx = [TDidx; TDidx_iti];
        %     end
        
        t3 = t2 + Duration;
        
        
        % Collect spikes/FR/rms for this stimulus/unit
        raster_x = [];
        raster_y = [];
        stim   = nan( numel(TDidx), Duration+1 );
        
        for it = 1:numel(TDidx)
                        
            stim(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) - 1;
            
            raster_x = [raster_x sp];
            raster_y = [raster_y it*ones(1,numel(sp))];
            
        end %it
        
        
        %% Add to plot
                keyboard % change to ticks
        subplot(hs(2)); hold on
        plot(raster_x, raster_y + add_y,...
            '.','MarkerSize',rasterdotsize,'Color',UnitColors(N_un+1,:))
        
        add_y = add_y + numel(TDidx);
        N_un = N_un + 1;
        
        
    end %iUn
    
    
    set(gca,'ylim',[0 add_y],'xlim',[0 Duration],'ytick',numel(TDidx))
    
    pd_x = 0:(1000/AMrates(stid-1)):Duration;
    pd_x = [pd_x; pd_x];
    pd_y = [zeros(1,size(pd_x,2)); add_y.*ones(1,size(pd_x,2))];
    
    plot(pd_x,pd_y,'k-')
    
    
    % Plot stimulus RMS
    subplot(hs(1)); hold on
    pd_y = [zeros(1,size(pd_x,2)); ones(1,size(pd_x,2))];
    plot(pd_x,pd_y,'k-')
    plot(0:Duration, mean(stim,1,'omitnan'),...
        'LineWidth',6,'Color',colors(stid,:))
    set(gca,'Color','none','xtick',[],'ytick',[])
    
    suptitle(Info.stim_ID_key{stid})
    
    savedir = '/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/ProcessedData/Rasters_selectUnPhases';
    savename = sprintf('SelectUnits_%s_%s_%s',SUBJECT,SESSION,Info.stim_ID_key{stid});
    print_eps_kp(hf(stid-1),fullfile(savedir,savename))
    
end %stid








end