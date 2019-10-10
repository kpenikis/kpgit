function MPH = makeMPHtable(TrialData,artifact_trs,spl,lpn,spiketimes,RateStream,OnlyGoodTransitions)
% MPH = makeMPHtable(TrialData,artifact_trs,spl,lpn,spiketimes,OnlyGoodTransitions)
% 
%   called by plotMPH_aggregate, MPHcontextComparisons 
% 
%   Does not limit/randomize number of trials. Do that later if desired.
% 
%   KP, 2019-04
% 

global AMrates rateVec_AC rateVec_DB trMin 


if nargin<7
    OnlyGoodTransitions=0;
end


% 
[Stream_FRsmooth,Stream_zscore] = convertSpiketimesToFR(round(spiketimes),TrialData.offset(end)+1,TrialData.onset(1),TrialData.offset(1),20,20,'silence');

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
all_TDidx = get_clean_trials(TrialData,artifact_trs,spl,lpn,0);
allStim = unique(TrialData.trID(all_TDidx));


%==========================================================
% Set up empty data matrices for this unit
%==========================================================

MPH = table;

MPH.AMrate       = 0;
MPH.PrevPd       = nan;
MPH.ThisStimID   = nan;
MPH.PrevStimID   = nan;
MPH.SeqPos       = nan;
MPH.Starttime    = nan;
MPH.Prev500msFR  = nan;
MPH.Prev100msFR  = nan;
MPH.PrevAMrt500  = nan;
MPH.PrevAMrt100  = nan;
MPH.raster       = {magic(3)};
MPH.x            = nan;
MPH.y            = nan;
MPH.FRsmooth     = {magic(3)};
MPH.zFR          = {magic(3)};
MPH.SPL          = 0;
MPH.LPN          = 0;

%==========================================================


%% Collect this unit's data

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% For each STIMULUS   (IR first)
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

for istim = allStim'
    
    if istim==1, continue, end   %no MPH for unmodulated
    
    
    %% 1) COLLECT AND SET STIMULUS INFO
    
    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
    
    
    % First get non-ITI trials
    
    st_TDidx = st_TDidx_ALL( TrialData.ITIflag(st_TDidx_ALL) == 0 );
    
    pst_TDidx = nan(size(st_TDidx));
    skip_it = [];
    
    if OnlyGoodTransitions
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
    end
    
    st_TDidx(skip_it)  = [];
    pst_TDidx = st_TDidx-1;
    
    
    % Add trials from ITI blocks
%     if istim==3 || istim==6
%         iti_TDidx = st_TDidx_ALL( TrialData.ITIflag(st_TDidx_ALL) == 1 );
%         skip_it = TrialData.trID(iti_TDidx-1)<6; %only use ITI trials that followed IR trials (thus, skip those after Warns)
%         iti_TDidx(skip_it) = [];
%         st_TDidx = [st_TDidx; iti_TDidx];
%         pst_TDidx = [pst_TDidx; iti_TDidx-1];
%     end
    
    
    % Skip transitions with too few trials
    prevStimIDs = unique(TrialData(pst_TDidx,:).trID)';
    transN = nan(1,numel(prevStimIDs));
    skip_idx=[];
    for pstid = prevStimIDs
        if numel(find(TrialData(pst_TDidx,:).trID==pstid)) < trMin
            skip_idx = [skip_idx find(pstid==prevStimIDs)];
        else
            transN(pstid==prevStimIDs) = numel(find(TrialData(pst_TDidx,:).trID==pstid));
        end
    end
    prevStimIDs(skip_idx)=[];
    transN(skip_idx)=[];
    
    if any(transN<trMin)
        keyboard
    end
    
    
    %% Set starttimes and period vectors
    
    clear t2
    t2 = TrialData.onset(st_TDidx);
    
    % Set the periods to gather
    switch istim
        case 7
            allPds = rateVec_AC;
        case 8
            allPds = rateVec_DB;
        otherwise
            allPds = repmat(AMrates(istim-1),1,ceil(AMrates(istim-1)));
    end
    
    
    
    %% 2) GET DATA
    
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . .
    % For each PREVIOUS stimulus
    % . . . . . . . . . . . . . . . . . . . . . . . . . . .
    
    for pstid = prevStimIDs
        
        % Create empty struct for gathering data
        MPH_temp = struct();
        for ipd = 1:numel(allPds)
            MPH_temp(ipd).x        = [];
            MPH_temp(ipd).y        = [];
            MPH_temp(ipd).raster   = zeros(transN(pstid==prevStimIDs),ceil(1000/allPds(ipd)));
            MPH_temp(ipd).FRsmooth = zeros(transN(pstid==prevStimIDs),ceil(1000/allPds(ipd)));
        end
        
        % Get TrialData indices for this transition 
        trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
        
        % Collect spikes/FR/rms for this transition
        for ipd = 1:numel(allPds)
            
            kt=0;
            for it = trans_TDidx'
                
                newPd_ts  =  t2(it)+[0 cumsum(1000./allPds)];
                if numel(newPd_ts)~=(numel(allPds)+1), keyboard, end
                
                
                % skip if this period extends beyond recording duration
                if ceil(newPd_ts(ipd))+ceil(1000/allPds(ipd))-1 > length(Stream_FRsmooth)
                    continue
                end
                
                kt=kt+1;
                                
                % Determine whether this period was in
                % the first or second half of the trial
                seqPos = (ipd>(length(allPds)/2))+1;
                if (seqPos<1 || seqPos>2)
                    keyboard
                end
                if (newPd_ts(ipd)-t2(it))<1 && ipd>1
                    keyboard
                end
                
                % Save info about this period
                MPH_temp(ipd).thisrate(kt)      = allPds(ipd);
                MPH_temp(ipd).seqPos            = seqPos;
                MPH_temp(ipd).pdtime(kt)        = newPd_ts(ipd)-t2(it);
                MPH_temp(ipd).Prev500msFR(kt)   = sum( spiketimes>=(newPd_ts(ipd)-500) & spiketimes<=newPd_ts(ipd))/500*1000;
                MPH_temp(ipd).Prev100msFR(kt)   = sum( spiketimes>=(newPd_ts(ipd)-100) & spiketimes<=newPd_ts(ipd))/100*1000;
                MPH_temp(ipd).PrevAMrt500       = getPrecRates(newPd_ts(ipd),500,RateStream);
                MPH_temp(ipd).PrevAMrt100       = getPrecRates(newPd_ts(ipd),100,RateStream);
                if ipd>1
                    MPH_temp(ipd).PrevPd        = allPds(ipd-1);
                else
                    prevTrID = TrialData(find([TrialData.onset]==t2(it))-1,:).trID;
                    if prevTrID~= pstid, keyboard, end
                    
                    if prevTrID < 7 && prevTrID >1 %if previous stim was pdc
                        MPH_temp(ipd).PrevPd    = AMrates(prevTrID-1);
                    elseif prevTrID ==7
                        MPH_temp(ipd).PrevPd    = rateVec_AC(end);
                    elseif prevTrID ==8
                        MPH_temp(ipd).PrevPd    = rateVec_DB(end);
                    end
                end
                
                % Put spikes into corresponding raster
                sp=[]; sp = spiketimes( spiketimes>(newPd_ts(ipd)) & spiketimes<=(newPd_ts(ipd+1)-1) ) - (newPd_ts(ipd));
                
                MPH_temp(ipd).raster(kt,ceil(sp))  = 1;
                MPH_temp(ipd).x                    = [MPH_temp(ipd).x sp];
                MPH_temp(ipd).y                    = [MPH_temp(ipd).y kt.*ones(size(sp))];
                MPH_temp(ipd).FRsmooth(kt,:)       = Stream_FRsmooth(ceil(newPd_ts(ipd))+[1:ceil(1000/allPds(ipd))]-1);
                MPH_temp(ipd).zFR(kt,:)            = Stream_zscore(ceil(newPd_ts(ipd))+[1:ceil(1000/allPds(ipd))]-1);
                
            end %it
            
            MPH_temp(ipd).raster((kt+1):end,:)   = [];
            MPH_temp(ipd).FRsmooth((kt+1):end,:) = [];
            
            if size(MPH_temp(ipd).raster,1) ~= size(MPH_temp(ipd).Prev500msFR,2)
                keyboard
            end
        end %ipd
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Put the data just collected into the table
        %    [ ThisStimID PrevStimID AMrate SeqPos Starttime Prev500msFR Prev100msFR (Prev500msPd) PrevPd raster spl lpn ]
        for ipd = 1:numel(allPds)
            
            if numel(unique(MPH_temp(ipd).PrevPd))>1 ...
                    || numel(unique(MPH_temp(ipd).seqPos))>1 ...
                    || numel(unique(round(MPH_temp(ipd).pdtime)))>1
                keyboard
            end
            
            MPH_addrow = { ...
                allPds(ipd)...
                mode(MPH_temp(ipd).PrevPd)  ...
                istim ...
                pstid  ...
                mode(MPH_temp(ipd).seqPos)  ...
                mode(round(MPH_temp(ipd).pdtime))...
                {MPH_temp(ipd).Prev500msFR} ...
                {MPH_temp(ipd).Prev100msFR} ...
                MPH_temp(ipd).PrevAMrt500 ...
                MPH_temp(ipd).PrevAMrt100 ...
                {MPH_temp(ipd).raster} {MPH_temp(ipd).x} {MPH_temp(ipd).y}...
                {MPH_temp(ipd).FRsmooth}  {MPH_temp(ipd).zFR}  spl  lpn };
            MPH = [MPH; MPH_addrow];
            
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    end %pstid (prev stim id)
    
    
    
end %istim


MPH(1,:) = [];   % mnNtrials = min(cell2mat(cellfun(@(x) size(x,1),MPH.raster,'UniformOutput',false)));



end

