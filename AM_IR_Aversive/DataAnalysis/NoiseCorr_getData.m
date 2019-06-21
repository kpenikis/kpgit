function NoiseCorr_getData
% 
%
% Would like to also calculate running correlation throughout session.
% See playNoiseCorr for example.
% But I don't think this method quantifies the right thing, what I'm after.
% What I want to know is if there are certain stimulus features that 
% induce a move toward higher or lower noise correlations. Not sure how to
% quantify that yet.
%
% KP, 2019-04
% 


%% Load Unit data files

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

subjcol     = [1 1 1];


%% Prepare figures
% 
% set(0,'DefaultTextInterpreter','none')
% set(0,'DefaultAxesFontSize',14)
% rng('shuffle')
% 
% scrsz = get(0,'ScreenSize');    %[left bottom width height]
% largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];
% tallrect  = [1 scrsz(4) scrsz(3)/3 scrsz(4)];
% 
% figwidthscales = [1.5 1 1 1 1 1 1.937 1.937];
% 
% 
% % Set colors
% colors = [   0 200 150;...
%     84  24  69;...
%     120  10  41;...
%     181   0  52;...
%     255  87  51;...
%     255 153   0]./255;
% colors = [ colors; ...
%     [37  84 156]./255 ;...
%     [19 125 124]./255 ];
% 
% rasterdotsize  = 18;
% psthlinewidth  = 4;


%% Step through Units

for iUn1 = 1:numel(UnitData)
        
    %%% skips merged units for now
    if numel(UnitInfo(iUn1,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
%     if UnitData(iUn).kw_p>0.05
%         continue
%     end
    
    subject     = UnitData(iUn1).Subject;
    session     = UnitData(iUn1).Session;
    channel     = UnitData(iUn1).Channel(1);
    Clu1        = UnitData(iUn1).Clu(1);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn1).spl;
    LP          = UnitData(iUn1).lpn;
    
    % Find other clus in this session
    theseUns = find(strcmp(subject,{UnitData.Subject}) & strcmp(session,{UnitData.Session}));
    theseUns = theseUns(theseUns>iUn1);
    
    for itU = 1:numel(theseUns)
        
        Clu2        = UnitData(theseUns(itU)).Clu(1);
        iUn2        = theseUns(itU);
        
        % Load data files
        if (iUn1>1 && ~( strcmp(subject,UnitData(iUn1-1).Subject) && strcmp(session,UnitData(iUn1-1).Session) )) || iUn1==1 && itU==1
            fprintf('Loading %s sess %s...\n',subject,session)
            clear TrialData Info RateStream SoundStream SpoutStream
            filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
            filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
        end
        if (iUn1>1 && ~( strcmp(subject,UnitData(iUn1-1).Subject) && strcmp(session,UnitData(iUn1-1).Session) && channel==UnitData(iUn1-1).Channel ) )  || iUn1==1  && itU==1
            clear Clusters Spikes
            filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
        end
        if ~exist('RateStream','var')
            keyboard
        end
        
        
        % Get spiketimes and shift based on calculated integration time
        if exist('Spikes','var')                                 % >>> UMS <<<
            
            spiketimes = unique(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==Clu1') * 1000 - spkshift);  %ms
            
        elseif exist('Clusters','var')                            % >>> KS <<<
            
            iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == Clu1);
            spiketimes = unique(Clusters(iClu).spikeTimes * 1000 - spkshift)';
            
        end
        
        fprintf(' analyzing clu %i X clu %i\n',Clu1,Clu2)
        
        
        %% Get raster data
        
%         convWin = 1;
        recLen = length(SoundStream);
        
        spiketimes1 = unique(round( Clusters([Clusters.clusterID]==Clu1).spikeTimes' * 1000 - spkshift ));
        sp1 = zeros(1,recLen);
        sp1(spiketimes1) = 1;
%         sp1 = convolveGauss(sp1,convWin);
        
        spiketimes2 = unique(round( Clusters([Clusters.clusterID]==Clu2).spikeTimes' * 1000 - spkshift ));
        sp2 = zeros(1,recLen);
        sp2(spiketimes2) = 1;
%         sp2 = convolveGauss(sp2,convWin);
        
        
        %% Now collect R noise values during each stimulus
        
        % Get all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials] = get_clean_trials( TrialData,...
            unique([Info.artifact(UnitData(iUn1).Channel).trials; Info.artifact(UnitData(iUn2).Channel).trials]),...
            UnitData(iUn1).spl,UnitData(iUn1).lpn);
        allStim = unique(TrialData.trID(all_TDidx))';
        
        minTrs = min(Ntrials(~isnan(Ntrials)));
        
        ymaxval = 0;
        
        Rnoise_RUN_stim     = nan(1,numel(allStim));
        Rnoise_RUN_stim_std = nan(1,numel(allStim));
        
        Rnoise_stim       = nan(1,numel(allStim));
        
        FR1_stim          = nan(1,numel(allStim));
        FR2_stim          = nan(1,numel(allStim));
        
        zFR1 = [];
        zFR2 = [];
        
        %% Collect trial indices and timestamps
        for ist = 1:numel(allStim)
            
            stid = allStim(ist);
            
            if stid==3 || stid==6
                TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
                % Find Pdc trials that follow same rate during ITI
                TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
                
                TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
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
                TDidx = [TDidx; TDidx_iti];
            end
            
            kt = 1:length(t2);
            t2     = t2(kt);
            TDidx  = TDidx(kt);
            
            t3 = t2 + Duration;
            
            
            %%
            % Preallocate
            Rstim  = nan( numel(TDidx), Duration+1 );
            FR1    = [];
            FR2    = [];
            
            % Collect rms/noise for this stimulus/unit
            for it = 1:numel(TDidx)
                
                Rstim(it,:) = ...
                    runningRcorr( (t2(it) : t3(it)) );
                
                FR1 = [FR1 sum(sp1((t2(it)+1):(t3(it))))];
                FR2 = [FR2 sum(sp2((t2(it)+1):(t3(it))))];
                
            end %it
            
            
            %% Store data
            % R noise -- raw trajectories
            Rstim;
            
            %  R noise -- norm
            Rstim' - repmat(Rstim(:,1)',size(Rstim,2),1);
            
            Rnoise_RUN_stim(ist)     = mean(mean(Rstim,2,'omitnan'),'omitnan');
            Rnoise_RUN_stim_std(ist) = std(mean(Rstim,2,'omitnan'),'omitnan');
            
            r = corrcoef(FR1,FR2);
            Rnoise_stim(ist)       = r(1,2);
            
            FR1_stim(ist) = mean(FR1);
            FR2_stim(ist) = mean(FR2);
            
            zFR1 = [zFR1 zscore(FR1)];
            zFR2 = [zFR2 zscore(FR2)];
            
        end %ist
        
        
        % Signal correlation
        r = corrcoef(FR1_stim,FR2_stim);
        plot([0 ist+1],[r(1,2) r(1,2)],'--k','LineWidth',2)
        
        % Noise correlation (all stim, method downer paper)
        r = corrcoef(zFR1,zFR2);
        plot([0 ist+1],[r(1,2) r(1,2)],'--m','LineWidth',2)
        
        % Noise correlation (whole session)
        minL = min([length(sp1) length(sp2)]);
        r = corrcoef(sp1(TrialData.onset(1):minL),sp2(TrialData.onset(1):minL));
        plot([0 ist+1],[r(1,2) r(1,2)],'--r','LineWidth',2)
        
        % Noise correlation (silence)
        r = corrcoef(sp1(TrialData.onset(1):TrialData.offset(1)),sp2(TrialData.onset(1):TrialData.offset(1)));
        plot([0 ist+1],[r(1,2) r(1,2)],'--b','LineWidth',2)
        
        % Each stimulus noise correlation - from running corr
        for ist = 1:numel(allStim)
            plot([allStim(ist); allStim(ist)],Rnoise_RUN_stim(ist)+Rnoise_RUN_stim_std(ist).*[-1; 1],...
                '-','Color',colors(allStim(ist),:),'LineWidth',4)
            plot(allStim(ist),Rnoise_RUN_stim(ist),'.','Color',colors(allStim(ist),:),'MarkerSize',50)
            plot(allStim(ist), Rnoise_stim(ist) ,'x','Color',colors(allStim(ist),:),'MarkerSize',20,'LineWidth',4)
        end
        
        
    end %iUn2
end %iUn1



end