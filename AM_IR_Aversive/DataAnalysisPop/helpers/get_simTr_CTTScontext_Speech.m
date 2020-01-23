function [SpikesTrials,StimTrials,includethiscell] = get_simTr_CTTScontext_Speech(UnitData,iUn)
% called by gatherCellTimeTrialStim_Speech
%
% KP, 2019-12-16
%

global fn trN Duration Stimuli spkshift k kfns SegDurs

includethiscell = 1;

spkshift = 0;

padPreT0 = 0;
padPost  = 0;

if ~exist('Stimuli','var')
    keyboard
    Stimuli = [1 5 6 2];
end
if ~exist('Duration','var')
    keyboard
    Duration  = 1000;
end

% Load data files
try
    fprintf('Loading %s sess %s...\n',UnitData(iUn).Subject,UnitData(iUn).Session)
    clear TrialData Clusters Spikes Info
    filename = sprintf( '%s_sess-%s_TrialData',UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Info',     UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,UnitData(iUn).Subject,UnitData(iUn).Session); load(fullfile(fn.processed,UnitData(iUn).Subject,filename));
    
    % Get spiketimes and shift based on calculated integration time
    if exist('Spikes','var')                                 % >>> UMS <<<
        spiketimes = unique(round(Spikes.sorted(UnitData(iUn).Channel).spiketimes(Spikes.sorted(UnitData(iUn).Channel).assigns==UnitData(iUn).Clu') * 1000 - spkshift));  %ms
    elseif exist('Clusters','var')                            % >>> KS <<<
        iClu = find([Clusters.maxChannel] == UnitData(iUn).Channel & [Clusters.clusterID] == UnitData(iUn).Clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
    end
catch
    keyboard
end


%%
% try

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials,~] = get_clean_trials(TrialData,unique(vertcat(Info.artifact(:).trials)),UnitData(iUn).spl,UnitData(iUn).lpn,1);
% [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(UnitData(iUn).Channel).trials,UnitData(iUn).spl,UnitData(iUn).lpn,1);

unStim = unique(TrialData.trID(all_TDidx))';

theseStim = intersect(Stimuli,unStim,'stable');

SpikesTrials  = nan(1,Duration+padPreT0+padPost,trN,numel(Stimuli),3);
StimTrials    = nan(1,Duration+padPreT0+padPost,trN,numel(Stimuli),3);

for ii = 1:numel(theseStim)
    
    stid = theseStim(ii);
    ist = find(stid==Stimuli);
    
    SegTemplate = getfield(k,kfns{ist});
        
    %% Collect trial indices and timestamps
    
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
    
    kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    t3 = TrialData.offset(TDidx);
    
    nReps = round(mean(t3-t2)./length(SegTemplate));
    
    tEdgS = repmat([0 length(SegTemplate)-1],nReps,1) + ([1:nReps]'-1)*(length(SegTemplate));
    
    
    %% Get spiking activity for each trial
    
    RMStr1st = nan(numel(t2),length(SegTemplate));
    RMStrRep = nan(numel(t2)*(nReps-1),length(SegTemplate));
    Rastr1st = zeros(numel(t2),length(SegTemplate));
    RastrRep = zeros(numel(t2)*(nReps-1),length(SegTemplate));
    
    iseg = 0;
    for it = 1:numel(t2)
        
        % First period
        sp = spiketimes( spiketimes>(t2(it) + tEdgS(1,1)) & spiketimes<=(t2(it) + tEdgS(1,2)) )  - (t2(it) + tEdgS(1,1));
        Rastr1st(it,sp) = 1;
        
        RMStr1st(it,:) = SoundStream( (t2(it) + tEdgS(1,1)) : (t2(it) + tEdgS(1,2)) );
        
        
        % Repeated periods
        if it==500
        figure;
        end
        for irep = 2:nReps
            iseg = iseg+1;
            
            % Find max xcorr with template (shift error accumulates if use
            % tEdgSeg directly
            leftEdge = t2(it) + tEdgS(irep,1)-10;
            nlags = 50;
            SearchRange = SoundStream( leftEdge + (0:length(SegTemplate)+19) );
            [c,lags] = xcov(SegTemplate,SearchRange/max(SearchRange),nlags);
            [~,im] = max(c);
            s1 = leftEdge-lags(im);
            
            if it==500
            clf
            plot(SegTemplate,'k','LineWidth',2)
            hold on
            plot( SoundStream(s1 + (0:length(SegTemplate)-1))/max(SoundStream(s1 + (0:length(SegTemplate)-1))) )
            pause(1)
            end
            
            RMStrRep(iseg,:) = SoundStream(s1 + (0:length(SegTemplate)-1));
            
            sp = spiketimes( spiketimes>s1 & spiketimes<=(s1 + length(SegTemplate)-1) )  - s1;
            RastrRep(iseg,sp) = 1; 
            
        end
        
    end %it
    
%     figure; hold on
%     plot(RMStrRep','b')
%     plot(RMStr1st','k')
    
    % Save data
    SpikesTrials(1,1:SegDurs(ist),1:size(Rastr1st,1),ist,1) = Rastr1st';
    SpikesTrials(1,1:SegDurs(ist),1:size(RastrRep,1),ist,2) = RastrRep';
    
    StimTrials (1,1:SegDurs(ist),1:1:size(RMStr1st,1),ist,1)  = RMStr1st';
    StimTrials (1,1:SegDurs(ist),1:1:size(RMStrRep,1),ist,2)  = RMStrRep';
    
    
    
    %% Now find this segment in the FULL SENTENCE context 
    
    stid_CONTEXT = [];
    switch stid
        case 1 %AsYou
            stid_CONTEXT = 4;
            appxT        = 4970;
        case 5 %Please
            stid_CONTEXT = 4;
            appxT        = 5360;
        case 6 %Trees
            stid_CONTEXT = 4;
            appxT        = 2170;
        case 2 %-ber
            stid_CONTEXT = 3;
            appxT        = 1700;
    end
    
    % Collect trial indices and timestamps
    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid_CONTEXT );
    
    % Get timestamps of onsets and offsets
    clear t2 t3 t_win
    t2 = TrialData.onset(TDidx);
    
    kt     = 1:length(t2);
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3 = TrialData.offset(TDidx);
    
    % Get spiking activity for each trial
    RMStrSen = nan(numel(t2),length(SegTemplate));
    RastrSen = zeros(numel(t2),length(SegTemplate));
    appxTime = nan(numel(t2),1);
    for it = 1:numel(t2)
        
        % Find the segment
        SearchRange = SoundStream( t2(it)+appxT + (0:length(SegTemplate)+50) );
        [c,lags] = xcov(SegTemplate,SearchRange/max(SearchRange));
        [~,im] = max(c);
        s1 = t2(it)+appxT-lags(im);
        appxTime(it) = s1-t2(it);
        
        RMStrSen(it,:) = SoundStream(s1 + (0:length(SegTemplate)-1));
        
        sp = spiketimes( spiketimes>s1 & spiketimes<=(s1 + length(SegTemplate)-1) )  - s1;
        RastrSen(it,sp) = 1;
        
    end %it
    
%     appxTime
%     plot(RMStrSen')
    
    % Save data
    SpikesTrials(1,1:SegDurs(ist),1:size(RastrSen,1),ist,3)   = RastrSen';
    StimTrials (1,1:SegDurs(ist),1:1:size(RMStrSen,1),ist,3)  = RMStrSen';
    
    
end %ist


end

