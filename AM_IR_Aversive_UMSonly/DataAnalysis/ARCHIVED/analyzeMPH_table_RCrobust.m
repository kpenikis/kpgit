function analyzeMPH_table_RCrobust(PlotThis)
%
%  MPHanalyses( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   Uses the TrialData (newer) version of saving stimulus info.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%   Option to plot only robust responders.
%
%  KP, 2018-04
%

% MPH consistency, reliability, etc 
%  as a function of previous stimulus
%  as a function of context 



global fn 
% close  all

if nargin<1
    PlotThis = 'FR'; 'FF'; 'TrV'; 'RS';
end

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
RobustUn =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  2;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
spktimeSHIFT = -0;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!



%% Set up figure


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
narrow = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

switch PlotThis
    case 'FR'
        yminval = 0;
        ymaxval = 40;
    case 'FF'
        yminval = 0;
        ymaxval = 4;
    case 'TrV'
        yminval = 0;
        ymaxval = 0.5;
    case 'RS'
        yminval = 0;
        ymaxval = 50;
    otherwise
        keyboard
end

AMrates = [2 4 8 16 32];
ContextStr = {'Periodic' 'Irregular'};
SeqPosStr = {'Early IR' 'Late IR'};

% Set up fig for population scatterplot
hf = figure;
set(hf,'Position',narrow)
for ir = 1:numel(AMrates)
for isp=1:2
    hsp(ir,isp)=subplot(numel(AMrates),2,ir*2-2+isp); hold on
    plot([yminval ymaxval],[yminval ymaxval],':k'); 
    axis square
    if (ir*2-2+isp)<3
        title(SeqPosStr{isp})
    end
    if (ir*2-2+isp)==9
        ylabel('Irregular context')
        xlabel('Periodic context')
    end
end
end
linkaxes(hsp,'xy')
xlim([yminval ymaxval])
ylim([yminval ymaxval])
suptitle([PlotThis ' for each period'])


colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        

%% Load Resp data table (with RCorr results)

fn = set_paths_directories('','',1);
Resp = readtable(fullfile(fn.processed,'RespTable_allSU_Rcorr'));

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;


N=0;

% Preallocate struct for results
foo = nan(0,1);
BinarySpikeData = cell(8,8);
[BinarySpikeData{:}] = deal(foo);

MTF_ts = [100 250 500 nan]';
foo = nan(8,8,3);

FR_Population=struct;
for ii=1:numel(MTF_ts)
    if isnan(MTF_ts(ii))
        MTF_ts_str{ii} = 'durFull';
    else
        MTF_ts_str{ii} = sprintf('dur%i',MTF_ts(ii));
    end
    FR_Population.(MTF_ts_str{ii}) = foo;
end


%%  SUBJECTS

% if nargin>0 && exist('select_subject','var')
%     subjects = {select_subject};
% else
    subjects = {'WWWf_253400' 'WWWlf_253395'};
% end

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

% if nargin>1 && exist('select_session','var')
%     Sessions = {select_session};
% else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
    Sessions = [];
    for ifn = 1:numel(SpkFns)
        if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
            Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
        end
    end
    
% end
% Sessions = flipud(Sessions);

% Step through each session
for sess = Sessions'
    
session = char(sess);


if RobustUn %skip sessions with no Robust SUs
    SessResp = Resp(strcmp(Resp.Subject,subject)&strcmp(Resp.Session,session),:);
    SessResp = SessResp(SessResp.AudResp_rc==1,:);
    if isempty(SessResp),  continue,   end
    Channels = unique(SessResp.ch)';
else 
    Channels =  [1:7 9:16];
end


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
    spikes = Spikes.sorted(channel);
    if RobustUn==0
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3) %no valid clus for this channel
            continue
        else
            Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
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
        spiketimes = spiketimes+spktimeSHIFT;
        
        if isempty(spiketimes)
            continue
        end
        
        
%         fprintf('-----------------------------\n')
        
        
        
        
        
        %% ################     NOW THE FUN STARTS     ##################
        
        
        
        
        
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Collect FR over session
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        bs_hist = 1;
        bs_smth = 20;
        
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');
        
        
        % Check if unit has very low overall firing rate
        if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
            continue
        end
        
        N = N+1;
        

        %%
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
%                 fprintf(' analyzing ch %i clu %i\n',channel,clu)
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
                all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                
                allStim = unique(TrialData.trID(all_TDidx));
                
                
                %==========================================================
                % Set up empty data matrices for this unit
                %==========================================================
                
                MPH = table;
                
                MPH.ThisStimID   = nan;
                MPH.PrevStimID   = nan;
                MPH.AMrate       = 0;
                MPH.SeqPos       = nan;
                MPH.Starttime    = nan;
                MPH.Prev500msFR  = nan;
%                 MPH.Prev500msPd  = nan;
%                 MPH.PrevPd       = nan;
                MPH.raster       = {magic(3)};
                MPH.SPL          = 0;
                MPH.LPN          = 0;
                
                %==========================================================
                
                
                %% Collect this unit's data
                
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                % For each STIMULUS   (IR first)
                % . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                                
                for istim = allStim'
                    
                    if istim==1, continue, end   %no MPH for unmodulated
                    
                    
                    %% 1) COLLECT AND SET SOME STIMULUS INFO 
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
                    
                    %%%  plot for Trial and ITI stimuli separately
                    ITIflag = unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for is = 1:numel(ITIflag)
                        
                        st_TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
                        
                        pst_TDidx = nan(size(st_TDidx));
                        skip_it = [];
                        
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
                        
                        % Now have final trials to work with
                        st_TDidx(skip_it)  = [];
                        pst_TDidx = st_TDidx-1;
                        
                        % Get duration/timestamps
                        clear t1 t2 t3 Durations t_win
                        t2 = TrialData.onset(st_TDidx);
                        t3 = TrialData.offset(st_TDidx);
                        Durations(2) = mode(diff([t2 t3],1,2));
                        t3 = t2 + Durations(2);
                        Durations(1) = mode(diff([(TrialData.onset(pst_TDidx)) t2],1,2));
                        t1 = t2 - Durations(1);
                        
                        % Set the periods to gather
                        switch istim
                            case 7
                                allPds = rateVec_AC;
                            case 8
                                allPds = rateVec_DB;
                            otherwise
                                allPds = repmat(AMrates(istim-1),1,ceil(AMrates(istim-1)/Durations(2)*1000));
                        end
                        
                        
                        
                        %% 2) GET DATA
                        
                        
%                         if istim<7
%                             % get corresponding periods
% %                             starttimes_to_match = unique(MPH(MPH.AMrate==AMrates(istim-1),:).Starttime)';
%                             these_rates = AMrates(istim-1);
%                         else 
%                             these_rates = AMrates;
%                         end  
                        
                        
                        % Get number of trials for each transition
                        transN = nan(1,numel(unique(TrialData(pst_TDidx,:).trID)));
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            transN(pstid==unique(TrialData(pst_TDidx,:).trID)') = numel(find(TrialData(pst_TDidx,:).trID==pstid));
                        end
                        
                        
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        % For each PREVIOUS stimulus
                        % . . . . . . . . . . . . . . . . . . . . . . . . . . .
                        
                        for pstid = unique(TrialData(pst_TDidx,:).trID)'
                            
                            
                            % Create empty struct for gathering data
                            MPH_temp = struct();
                            for ipd = 1:numel(allPds)
                                MPH_temp(ipd).raster  = zeros(transN(pstid==unique(TrialData(pst_TDidx,:).trID)'),ceil(1000/allPds(ipd)));
                            end
                            
                            % Get TrialData indices for this transition and shuffle
                            trans_TDidx = find(TrialData(pst_TDidx,:).trID==pstid);
                            trans_TDidx = trans_TDidx(randperm(numel(trans_TDidx)));
                            
                            % Collect spikes/FR/rms for this transition
                            for it = 1:numel(trans_TDidx)
                                
%                                 if it>min(transN), continue, end  %truncate N trials at min
                                	newPds = t2(it)+[0 cumsum(1000./allPds)];
                                    
                                for ipd = 1:numel(allPds)
                                    
                                    sp=[]; sp = spiketimes( spiketimes>round(newPds(ipd)) & spiketimes<=round(newPds(ipd+1)-1) ) - round(newPds(ipd));
                                    
                                    % Determine whether this period was in
                                    % the first or second half of the trial
                                    seqPos = (ipd>(length(allPds)/2))+1;
                                    if (seqPos<1 || seqPos>2)
                                        keyboard
                                    end
                                    
                                    
                                    % Save some info about this period
%                                     MPH_temp(these_rates==allPds(ipd)).seqPos(it,seqPos)  = seqPos;
%                                     MPH_temp(these_rates==allPds(ipd)).pdtime(it,seqPos)  = newPds(ipd)-t2(it);
%                                     MPH_temp(these_rates==allPds(ipd)).Prev500msFR(it,seqPos) = sum( spiketimes>=(newPds(ipd)-500) & spiketimes<=newPds(ipd))/500;
                                    MPH_temp(ipd).thisrate(it)= allPds(ipd);                                    
                                    MPH_temp(ipd).seqPos      = seqPos;
                                    MPH_temp(ipd).pdtime(it)  = newPds(ipd)-t2(it);
                                    MPH_temp(ipd).Prev500msFR(it) = sum( spiketimes>=(newPds(ipd)-500) & spiketimes<=newPds(ipd))/500*1000;
                                    
                                    % Put spikes into corresponding raster
%                                     MPH_temp(AMrates==thesePds(ipd)).raster(min(transN)*(seqPos-1)+it,sp) = 1;
                                    MPH_temp(ipd).raster(it,sp) = 1;
                                    
                                    
                                end %ipd
                                
                            end %it
                            
                            
                            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            % Put the data just collected into the table
                            %    [ AMrate ThisStimID PrevStimID SeqPos Starttime Prev500msFR   * Prev500msPd PrevPd *   raster spl lpn ]
                            for ipd = 1:numel(allPds)
                                    
                                    Prev500msFR = mean(MPH_temp(ipd).Prev500msFR);
%                                     Prev500msPd = [];
%                                     PrevPd
                                    seqPos    = mode(MPH_temp(ipd).seqPos);
                                    starttime = mode(round(MPH_temp(ipd).pdtime));
                                    
                                    MPH_addrow = {istim pstid allPds(ipd) seqPos starttime Prev500msFR  {MPH_temp(ipd).raster} spl lpn};
                                    MPH = [MPH; MPH_addrow];
                                    
                            end
                            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            
                            clear MPH_temp                            
                            
                            
                        end %pstid (prev stim id)
                        
                        
                    end %ITIflag (Trial or ITI)                    
                end %istim
                
                MPH(1,:) = [];
                
                
                %%  3) Analyze data
                
                % Separate periods by AM rate, prevStim irrelevant
                % Looking for difference between contexts in:
                % (measure within context, for each prevStim then avg'd)
                % IR * seqPos, matched from Pdc
                %   FR
                %   FF
                %   variability of trial-trial responses (binned?)
                %   VS
                % (measure across contexts)
                %   ** correlation
                %
                
                
                % Later, this will be selected from this unit's BMF or
                % something from Resp table
                anThisRate = [2 4 8 16 32];
                
                for this_rate = anThisRate
                    
                    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % Get IR context first, and collect starttimes to match
                    
                    Period = 1000/this_rate;
                    theseIRs = unique(MPH.ThisStimID(MPH.ThisStimID>6))';
                    PdData = struct;
                    
                    for IRstim = theseIRs
                        
                        subMPH = MPH(MPH.ThisStimID==IRstim & MPH.AMrate==this_rate,:);
                        
                        if isempty(subMPH), continue, end
                                                
                        for iseq = 1:2
                            
                            % Preallocate
                            FR  = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            SEM = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            FF  = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            TrV = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            VS  = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            RS  = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            RP  = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            
                            for ipst = 1:numel(subMPH(subMPH.SeqPos==iseq,:).raster)
                                
                                raster = subMPH(subMPH.SeqPos==iseq,:).raster{ipst};
                                if ~(ceil(Period)==size(raster,2))
                                    keyboard
                                end
                                
                                FR(ipst)  = sum(sum( raster ))/size(raster,1)/Period*1000;
                                SEM(ipst) = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                                FF(ipst)  = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                                TrV(ipst) = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                                
                                % Reformat raster for VS calculation
                                spktimes=[];
                                for it = 1:size(raster,1)
                                    spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                                end
                                [VS(ipst),RS(ipst),RP(ipst),] = vectorstrength(spktimes,Period);
                                
                            end
                            
                            % Average over previous stimuli and save for
                            % comparison
                            
                            idx = 2*(IRstim-7) + iseq;
                            
                            PdData(idx,1).stimID = IRstim;
                            PdData(idx,1).seqPos = iseq;
                            PdData(idx,1).starttime = mode(subMPH(subMPH.SeqPos==iseq,:).Starttime);
                            PdData(idx,1).FR  = mean(FR,'omitnan');
                            PdData(idx,1).SEM = mean(SEM,'omitnan');
                            PdData(idx,1).FF  = mean(FF,'omitnan');
                            PdData(idx,1).TrV = mean(TrV,'omitnan');
                            PdData(idx,1).VS  = mean(VS,'omitnan');
                            PdData(idx,1).RS  = mean(RS,'omitnan');
                            PdData(idx,1).RP  = mean(RP,'omitnan');
                            
                            
                        end %iseq
                        
                        
                    end %istim -- IRREGULAR first
                    
                    
                    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % Now get matching PERIODIC context data
                    
                    subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate,:);
                    pdcstarttimes = unique(subMPH.Starttime)';
                    
                    for iIRpd = 1:size(PdData,1)
                        
                        if isempty(PdData(iIRpd).stimID)
                            continue
                        end
                        
                        IRstim = PdData(iIRpd).stimID;
                        iseq = PdData(iIRpd).seqPos;
                        idx = 2*(IRstim-7) + iseq;
                        if idx~=iIRpd, keyboard, end
                        
                        [~,ipd] = min(abs(PdData(idx).starttime - pdcstarttimes));
                        
                        % Preallocate
                        FR  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        SEM = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        FF  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        TrV = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        VS  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        RS  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        RP  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        
                        for ipst = 1:numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster)
                            
                            raster = subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster{ipst};
                            if ~(ceil(Period)==size(raster,2))
                                keyboard
                            end
                            
                            FR(ipst)   = sum(sum( raster ))/size(raster,1)/Period*1000;
                            SEM(ipst)  = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                            FF(ipst)   = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                            TrV(ipst)  = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                            
                            % Reformat raster for VS calculation
                            spktimes=[];
                            for it = 1:size(raster,1)
                                spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                            end
                            [VS(ipst),RS(ipst),RP(ipst),] = vectorstrength(spktimes,Period);
                            
                        end %ipst
                        
                        
                        PdData(idx,2).stimID = this_rate;
                        PdData(idx,2).seqPos = iseq;
                        PdData(idx,2).starttime = mode(subMPH(subMPH.SeqPos==iseq,:).Starttime);
                        PdData(idx,2).FR  = mean(FR,'omitnan');
                        PdData(idx,2).SEM = mean(SEM,'omitnan');
                        PdData(idx,2).FF  = mean(FF,'omitnan');
                        PdData(idx,2).TrV = mean(TrV,'omitnan');
                        PdData(idx,2).VS  = mean(VS,'omitnan');
                        PdData(idx,2).RS  = mean(RS,'omitnan');
                        PdData(idx,2).RP  = mean(RP,'omitnan');
                        
                        
                        
                        %% And add these points to the population plot
                        
                        subplot(hsp(AMrates==this_rate,iseq)); hold on
                        
                        switch PlotThis
                            
                            case 'FR'                                
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,2).SEM,'horizontal','Color',colors(subMPH.ThisStimID(1),:))
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,1).SEM,'vertical','Color',colors(subMPH.ThisStimID(1),:))
                                plot(PdData(idx,2).FR, PdData(idx,1).FR,'o',...
                                    'Color',colors(subMPH.ThisStimID(1),:),'LineWidth',1.5,'MarkerSize',10)
                                
                            case 'FF'
                                plot(PdData(idx,2).FF, PdData(idx,1).FF,'o',...
                                    'Color',colors(subMPH.ThisStimID(1),:),'LineWidth',1.5,'MarkerSize',10)
                                
                            case 'TrV'
                                plot(PdData(idx,2).TrV, PdData(idx,1).TrV,'o',...
                                    'Color',colors(subMPH.ThisStimID(1),:),'LineWidth',1.5,'MarkerSize',10)
                                
                            case 'RS'
                                plot(PdData(idx,2).RS, PdData(idx,1).RS,'o',...
                                    'Color',colors(subMPH.ThisStimID(1),:),'LineWidth',1.5,'MarkerSize',10)  
                                
                        end
                        hold off
                        
                        
                        
                    end %iIRpd
                    
                    
                    
                    
%                     for ipd = 1:size(PdData,1)
%                         
%                         errorbar(PdData(ipd,2).FR, PdData(ipd,1).FR, PdData(ipd,2).SEM,'horizontal')
%                         errorbar(PdData(ipd,2).FR, PdData(ipd,1).FR, PdData(ipd,1).SEM,'vertical')
%                         plot(PdData(ipd,2).FR, PdData(ipd,1).FR,'o')
                        
% %                         figure;
% % %                         errorbar([PdData(:,2).seqPos]-0.1,[PdData(:,2).FR],[PdData(:,2).SEM],...
% % %                             '.k','LineWidth',2)
% %                         hold on
% %                         bar([PdData(:,2).seqPos]-0.1,[PdData(:,2).TrV],...
% %                             'EdgeColor','k','FaceColor','none','LineWidth',2,'BarWidth',0.5)
% % %                         errorbar([PdData(:,1).seqPos]+0.1,[PdData(:,1).FR],[PdData(:,1).SEM],...
% % %                             '.b','LineWidth',2)
% %                         bar([PdData(:,1).seqPos]+0.1,[PdData(:,1).TrV],...
% %                             'EdgeColor','b','FaceColor','none','LineWidth',2,'BarWidth',0.5)
% %                         xlim([0 3])
% %                         
% %                         return
% %                         
% %                         errorbar(PdData(ipd,2).seqPos,PdData(ipd,2).FR,PdData(ipd,2).SEM)
                        
%                         PdData(ipd,:).SEM
%                         PdData(ipd,:).FF
%                         PdData(ipd,:).VS
%                         PdData(ipd,:).RS
%                         PdData(ipd,:).RP
%                     end
                    
% %                     
                    
                end  %this_rate
                
                
                %%   End of Datapoint
                                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects


%% Save figure(s)

savedir = fullfile(fn.processed,'MPHmetrics');
if SUonly==1
    savedir = [savedir '/SU'];
end
if RobustUn==1
    savedir = [savedir '/RobustOnly'];
end
if ~exist(savedir,'dir')
    mkdir(savedir);
end
savename = sprintf('AllMPHs_Context_%s',PlotThis);
print_eps_kp(hf,fullfile(savedir,savename))








end %function




