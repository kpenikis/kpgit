function [UnitInfo, UnitData, Info, TrialData, Clusters, StimResp, artifactTrs ] = collectRasterDataSession(SUBJECT,SESSION)
%
% [UnitInfo, UnitData, Info, TrialData, Clusters, StimResp, artifactTrs ] = collectRasterDataSession(SUBJECT,SESSION)
%
%   ONLY ONE SESSION CURRENTLY
% 
%   Collects data for all stimuli and all units for further analyses.
%   Trials restricted to only those labeled clean for all relevant
%   channels.
%   Currently set to include ITI trials (for 4 and 32 Hz).
%   The resulting raster matrices are output in two places.
%   Organized by stimulus:
%       StimResp(stid).raster        (trials,spikes,Units)
%       StimResp(stid).trIdx
%   And organized by unit:
%       UnitData(iUn).raster{stid}   (trials,spikes)
%   The order of units is the same in both!
%
% KP, 2019-03
%

global RateStream

%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));


% Filter Unit files to just this session and sort by baseline FR
UnitData = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
[~, baseFRrank] = sort([UnitData.BaseFR]);
UnitInfo = UnitInfo(baseFRrank,:);
UnitData = UnitData(baseFRrank);

% Get indices of narrow spikes and regular spiking units
UnitInfo = labelRSNS(UnitInfo);

% Also find trials to skip
Channels = unique([UnitInfo.Channel]);
artifactTrs = [];
for ich = Channels'
    artifactTrs = [artifactTrs Info.artifact(ich).trials'];
end
artifactTrs = unique(artifactTrs);

minTrs = 12;


%% Find stimuli and trial indices to include in analysis

% Get sound parameters
[dBSPL,LP] = theseSoundParams(TrialData);
if numel(dBSPL)>1 || numel(LP)>1
    keyboard
end

% Find all stimuli presented with these parameters, given a sufficient minimum
% number of trials without disruptive artifact, while the animal was drinking
[all_TDidx,Ntrials,~,allStim] = get_clean_trials(TrialData,artifactTrs,dBSPL,LP);

if sum(Ntrials < minTrs)==1
    keyboard
    all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
    allStim(Ntrials<minTrs)  = [];
    Ntrials(Ntrials<minTrs) = [];
elseif  sum(Ntrials < minTrs)>1
    keyboard
end


% Now step through stimuli and units
StimResp = struct();

for stid = unique([TrialData.trID])'
    
    if stid==0, continue, end
    
    for iUn = 1:size(UnitInfo,1)
        
        %%% still must edit for merged units
        %     if numel(UnitInfo(iUn,:).Session{:})>2  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        %         continue
        %     end
        
        channel = UnitData(iUn).Channel(1);
        clu     = UnitData(iUn).Clu(1);
        
        % Get spiketimes (KS)
        spiketimes = round(Clusters(([Clusters.maxChannel]==channel & [Clusters.clusterID]==clu)).spikeTimes*1000 - spkshift)';
        
        [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,50,'silence');
        
        
        %%
        % Get trial numbers for this stimulus
        
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
        
        
        % Timestamps of onsets and offsets
        
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        Duration = mode(diff([t2 t3],1,2));
        
        % Add ITI trials (shortened to match duration)
        if ~isempty(TDidx_iti)
            t2 = [t2; TrialData.onset(TDidx_iti)];
            TDidx = [TDidx; TDidx_iti];
        end
        
        t3 = t2 + Duration;
        
        
        % Collect spikes/FR for this stimulus/unit
        rastermat = zeros( numel(TDidx), Duration );
        for it = 1:numel(TDidx)
            
            sp=[]; sp = spiketimes( spiketimes>t2(it) ...
                & spiketimes<=t3(it) ) - t2(it);
            
            rastermat(it,sp) = 1;
            
        end %it
        
        
        %% Store unit raster
        
        UnitData(iUn).raster{stid} = rastermat;
        StimResp(stid).raster(:,:,iUn) = rastermat;
        StimResp(stid).trIdx = TDidx;
        
        
    end %iUn
    
end %stid

end
