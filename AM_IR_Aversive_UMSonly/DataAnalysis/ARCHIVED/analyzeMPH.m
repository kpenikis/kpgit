function analyzeMPH(PlotThis,select_stim)
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
close  all

if nargin<1
    PlotThis = 'FR'; 'FF'; 'TrV'; 'VS'; 'MP'; 'CV'; 'r_all'; 'r_pdc';
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
elseif nargin==1
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
end

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
AMrespFilt =  1;
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
narrow     = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
smallrect  = [1 scrsz(4)/3 scrsz(3)/5 scrsz(4)/3];
largerrect = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];

switch PlotThis
    case 'FR'
        yminval = 0;
        ymaxval = 40;
    case 'FF'
        yminval = 0;
        ymaxval = 10;
    case 'TrV'
        yminval = 0;
        ymaxval = 0.7;
    case 'CV'
        yminval = 0;
        ymaxval = 4;
    case 'VS'
        yminval = 0;
        ymaxval = 1;
    case 'MP'
        yminval = 0;
        ymaxval = 2*pi;
    case {'r_all' 'r_pdc' 'avTrCorr'}
        yminval = 0;
        ymaxval = 1;
    otherwise
        keyboard
end

AMrates = [2 4 8 16 32];
ContextStr = {'Periodic' 'Irregular'};
SeqPosStr = {'Early' 'Late'};

% Set up figs for population scatterplots

% Each rate and time separately
hfe = figure;
set(hfe,'Position',narrow)
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
hold off


% Each rate separately, time merged
hfr = figure;
set(hfr,'Position',narrow)
for ir = 1:numel(AMrates)
    hspr(ir)=subplot(numel(AMrates),1,ir); hold on
    plot([yminval ymaxval],[yminval ymaxval],':k'); 
    axis square
    title([num2str(AMrates(ir)) ' Hz'])
    if ir==5
        ylabel('Irregular context')
        xlabel('Periodic context')
    end
end
linkaxes(hspr,'xy')
xlim([yminval ymaxval])
ylim([yminval ymaxval])
hold off


% All rates together
hfo = figure; hold on
set(hfo,'Position',smallrect)
axis square
plot([yminval ymaxval],[yminval ymaxval],':k');
xlim([yminval ymaxval])
ylim([yminval ymaxval])
ylabel('Irregular context')
xlabel('Periodic context')
hold off


% All rates together - Early/Late
hfoel = figure; hold on
set(hfoel,'Position',largerrect)
for isp=1:2
    hsel(isp) = subplot(1,2,isp); hold on
    axis square
    plot([yminval ymaxval],[yminval ymaxval],':k');
    xlim([yminval ymaxval])
    ylim([yminval ymaxval])
    ylabel('Irregular context')
    xlabel('Periodic context')
    title(SeqPosStr{isp})
    hold off
end
hold off



if strcmp(select_stim,'all')
    
    % Diff IR-Pdc, as a function of AMrate
    hfdrData=[];
    hfdr = figure; hold on
    set(hfdr,'Position',narrow)
    plot([0 6],[0 0],':k');
    xlim([0 6])
    ylim((ymaxval)*[-1 1])
    ylabel('diff IR-Pdc')
    xlabel('AM rate')
    set(gca,'xtick',1:5,'xticklabel',{'2 Hz' '4 Hz' '8 Hz' '16 Hz' '32 Hz'})
    hold off
    
    % Diff IR-Pdc, as a function of distance to BMF-fr
    hfddData=[];
    hfdd = figure; hold on
    set(hfdd,'Position',narrow)
    plot([-5 5],[0 0],':k');
    xlim([-5 5])
    ylim((ymaxval)*[-1 1])
    ylabel('diff IR-Pdc')
    xlabel('log distance to BMF')
    hold off
    
end

% % MTF
% hfmtf = figure; hold on
% set(hfmtf,'Position',fullscreen)
% xlim([0 6])
% ylim([yminval ymaxval])
% ylabel(PlotThis)
% xlabel('AM rate')
% set(gca,'xtick',1:5,'xticklabel',{'2 Hz' '4 Hz' '8 Hz' '16 Hz' '32 Hz'})
% hold off

% Regress IR-Pdc difference against previous period, prev FR
hfppData=[];




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
        
alphval = 0.6;
dotsize = 120;
        

%% Load Resp data table (with RCorr results)

fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'RespStruct_allSU'));
Resp = q.Resp;

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;


N=0;

% Preallocate struct for results
FinalDPs=struct;
for ir=1:5
    for is=1:2
        FinalDPs(ir,is).Pdc = [];
        FinalDPs(ir,is).IR  = [];
        DP_TrV(ir,is).Pdc = [];
        DP_TrV(ir,is).IR  = [];
        DP_FR(ir,is).Pdc  = [];
        DP_FR(ir,is).IR   = [];
    end
end


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


if AMrespFilt %skip sessions with no Robust SUs
    SessResp = Resp(strcmp({Resp.Subject},subject)&strcmp({Resp.Session},session));
    if isempty(SessResp),  continue,   end
    
    % Here add call to a function that outputs the indices of SessResp 
    % that contain a responsive unit. With this method, the thresholds can
    % be flexibly changed.
    [sig_ids,SessResp] = identifyResponsiveUnits(SessResp);
    
    Channels = unique([SessResp(sig_ids).Channel]);
    
    if isempty(Channels), continue, end
    
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
    if AMrespFilt==0
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3) %no valid clus for this channel
            continue
        else
            Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
        end
    else
        Clus = SessResp([SessResp.Channel]==channel).Clu;
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
        
        ThisResp = SessResp([SessResp.Channel]==channel & [SessResp.Clu]==clu);
        
        if mean(mean(ThisResp.FR_raw_tr,1,'omitnan'),2,'omitnan')>8
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
        
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(spiketimes,...
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
                MPH.Prev100msFR  = nan;
                MPH.PrevPd       = nan;
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
                    ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
                    
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
                                    MPH_temp(ipd).thisrate(it)= allPds(ipd);                                    
                                    MPH_temp(ipd).seqPos      = seqPos;
                                    MPH_temp(ipd).pdtime(it)  = newPds(ipd)-t2(it);
                                    MPH_temp(ipd).Prev500msFR(it) = sum( spiketimes>=(newPds(ipd)-500) & spiketimes<=newPds(ipd))/500*1000;
                                    MPH_temp(ipd).Prev100msFR(it) = sum( spiketimes>=(newPds(ipd)-100) & spiketimes<=newPds(ipd))/100*1000;
                                    if ipd>1
                                        MPH_temp(ipd).PrevPd       = allPds(ipd-1);
                                     else
                                        prevTrID = TrialData(find([TrialData.onset]==t2(it))-1,:).trID;
                                        if prevTrID < 7 && prevTrID >1 %if previous stim was pdc
                                            MPH_temp(ipd).PrevPd   = AMrates(prevTrID-1);
                                        elseif prevTrID ==7
                                            MPH_temp(ipd).PrevPd   = rateVec_AC(end);
                                        elseif prevTrID ==8
                                            MPH_temp(ipd).PrevPd   = rateVec_DB(end);
                                        end
                                    end
                                    
                                    
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
                                    Prev100msFR = mean(MPH_temp(ipd).Prev100msFR);
                                    PrevPd      = mode(MPH_temp(ipd).PrevPd);
                                    seqPos      = mode(MPH_temp(ipd).seqPos);
                                    starttime   = mode(round(MPH_temp(ipd).pdtime));
                                    
                                    MPH_addrow = {istim pstid allPds(ipd) seqPos starttime Prev500msFR Prev100msFR PrevPd  {MPH_temp(ipd).raster} spl lpn};
                                    MPH = [MPH; MPH_addrow];
                                    
                            end
                            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            
                            clear MPH_temp                            
                            
                            
                        end %pstid (prev stim id)
                        
                        
                    end %ITIflag (Trial or ITI)                    
                end %istim
                
                MPH(1,:) = [];
                
%                 mnNtrials = min(cell2mat(cellfun(@(x) size(x,1),MPH.raster,'UniformOutput',false)));
                
                


                
                %%  3) ANALYZE DATA
                
                
                % Rates to plot are selected by input variable select_stim,
                % using data from the Resp table. For example, can choose
                % to plot only this unit's BMF.
                
                if strcmp(select_stim,'all')
                    ThisResp.all = 1:5;
                    diffData_unit= [];
                end
                MTFdata = [];
                if ~isfield(ThisResp,select_stim)
                    continue
                end
                
                
                if strcmp(PlotThis,'VS') && strcmp(select_stim,'iBMF_FR')
                    if ~ismember(ThisResp.(select_stim), ThisResp.iSync)
%                         keyboard
                    end
                end
                
                
                anThisRate = AMrates(ThisResp.(select_stim));
                
                for this_rate = anThisRate
                                        
                    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % Get average periodic MPH
                    
                    MPHthisRate = MPH(MPH.AMrate==this_rate,:);
                    
                    iPdc = find(MPHthisRate.ThisStimID<7)';
                    iPdc  = randperm(length(iPdc));
                    iIR  = find(MPHthisRate.ThisStimID>=7)';
                    iIR  = randperm(length(iIR))+max(iPdc);
                    npds = min(length(iPdc),length(iIR));
                    
                    try
                    raster_AllPds = vertcat(MPHthisRate([iPdc(1:npds) iIR(1:npds)],:).raster{:});
                    MPH_AllPds = sum(raster_AllPds,1)./size(raster_AllPds,1)*1000;
                    
                    raster_PdcPds = vertcat(MPH(MPH.AMrate==this_rate & MPH.ThisStimID<7,:).raster{:});
                    MPH_PdcPds = sum(raster_PdcPds,1)./size(raster_PdcPds,1)*1000;
                    catch
                        keyboard
                    end
                    
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
                            FR       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            SEM      = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            FF       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            TrV      = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            CV       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            VS       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            VSp      = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            MP       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            rMPH_all = nan(2,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            rMPH_pdc = nan(2,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                            
                            for ipst = 1:numel(subMPH(subMPH.SeqPos==iseq,:).raster)
                                
                                raster = subMPH(subMPH.SeqPos==iseq,:).raster{ipst};
                                if ~(ceil(Period)==size(raster,2))
                                    keyboard
                                end
                                
                                FR(ipst)  = sum(sum( raster ))/size(raster,1)/Period*1000; 
                                SEM(ipst) = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                                FF(ipst)  = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                                TrV(ipst) = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                                CV(ipst)  = calc_binned_bootstrapped_CV(raster,tVarBin);
                                [r,p] = corrcoef(MPH_AllPds,sum(raster,1)/size(raster,1)*1000);
                                rMPH_all(:,ipst) = [r(1,2); p(1,2)];
                                [r,p] = corrcoef(MPH_PdcPds,sum(raster,1)/size(raster,1)*1000);
                                rMPH_pdc(:,ipst) = [r(1,2); p(1,2)];
                                
                                % Reformat raster for VS calculation
                                spktimes=[];
                                for it = 1:size(raster,1)
                                    spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                                end
                                [VS(ipst),~,VSp(ipst)] = vectorstrength(spktimes,Period);
                                
                                % Also calculate mean phase
                                [MP(ipst),R(ipst)] = meanphase(spktimes,Period);
                                
                                
                            end
                            
%                             % Remove VS if sync ns
%                             if mean(VSp,'omitnan') > 0.001
%                                 VS = 0;
%                             end
                            
                            
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
                            PdData(idx,1).CV  = mean(CV,'omitnan');
                            PdData(idx,1).VS  = mean(VS,'omitnan');
                            PdData(idx,1).MP  = mean(MP,'omitnan');
                            PdData(idx,1).r_all  = mean(rMPH_all(1,:));
                            PdData(idx,1).p_all  = mean(rMPH_all(2,:));
                            PdData(idx,1).r_pdc  = mean(rMPH_pdc(1,:));
                            PdData(idx,1).p_pdc  = mean(rMPH_pdc(2,:));
                            PdData(idx,1).prFR500 =  mean(subMPH(subMPH.SeqPos==iseq,:).Prev500msFR);
                            PdData(idx,1).prFR100 =  mean(subMPH(subMPH.SeqPos==iseq,:).Prev100msFR);
                            PdData(idx,1).prPd =  mode(subMPH(subMPH.SeqPos==iseq,:).PrevPd);
                                                        
                            
                        end %iseq
                        
                        
                    end %IRstim -- IRREGULAR first
                    
                    
                    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    % Now get matching PERIODIC context data
                    
                    subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate,:);
                    pdcstarttimes = unique(subMPH.Starttime)';
                    
                    for iIRpd = 1:size(PdData,1)
                        
                        if isempty(PdData(iIRpd,1).stimID)
                            continue
                        end
                        
                        IRstim = PdData(iIRpd,1).stimID;
                        iseq = PdData(iIRpd,1).seqPos;
                        idx = 2*(IRstim-7) + iseq;
                        if idx~=iIRpd, keyboard, end
                        
                        strtdiffs = PdData(idx,1).starttime - pdcstarttimes;
                        [~,ipd] = min(strtdiffs(strtdiffs>=0));
                        
                        % Preallocate
                        FR  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        SEM = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        FF  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        TrV = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        CV  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        VS  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        VSp = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        MP  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        rMPH_all = nan(2,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        rMPH_pdc = nan(2,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
                        
                        for ipst = 1:numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster)
                            
                            raster = subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster{ipst};
                            if ~(ceil(Period)==size(raster,2))
                                keyboard
                            end
                            
                            FR(ipst)   = sum(sum( raster ))/size(raster,1)/Period*1000;
                            SEM(ipst)  = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                            FF(ipst)   = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                            TrV(ipst)  = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                            CV(ipst)   = calc_binned_bootstrapped_CV(raster,tVarBin);
                            [r,p] = corrcoef(MPH_AllPds,sum(raster,1)/size(raster,1)*1000);
                            rMPH_all(:,ipst) = [r(1,2); p(1,2)];
                            [r,p] = corrcoef(MPH_PdcPds,sum(raster,1)/size(raster,1)*1000);
                            rMPH_pdc(:,ipst) = [r(1,2); p(1,2)];
                            
                            % Reformat raster for VS calculation
                            spktimes=[];
                            for it = 1:size(raster,1)
                                spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                            end
                            [VS(ipst),~,VSp(ipst)] = vectorstrength(spktimes,Period);
                            
                            % Also calculate mean phase
                            [MP(ipst),R(ipst)] = meanphase(spktimes,Period);
                            
                        end %ipst
                        
                        
%                         % Remove VS if sync ns
%                         if mean(VSp,'omitnan') > 0.001
%                             VS = 0;
%                         end
                        
                        
                        PdData(idx,2).stimID = this_rate;
                        PdData(idx,2).seqPos = iseq;
                        PdData(idx,2).starttime = mode(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Starttime);
                        PdData(idx,2).FR  = mean(FR,'omitnan');
                        PdData(idx,2).SEM = mean(SEM,'omitnan');
                        PdData(idx,2).FF  = mean(FF,'omitnan');
                        PdData(idx,2).TrV = mean(TrV,'omitnan');
                        PdData(idx,2).CV  = mean(CV,'omitnan');
                        PdData(idx,2).VS  = mean(VS,'omitnan');
                        PdData(idx,2).MP  = mean(MP,'omitnan');
                        PdData(idx,2).r_all  = mean(rMPH_all(1,:));
                        PdData(idx,2).p_all  = mean(rMPH_all(2,:));
                        PdData(idx,2).r_pdc  = mean(rMPH_pdc(1,:));
                        PdData(idx,2).p_pdc  = mean(rMPH_pdc(2,:));
                        PdData(idx,2).prFR500 =  mean(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Prev500msFR);
                        PdData(idx,2).prFR100 =  mean(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Prev100msFR);
                        PdData(idx,2).prPd    =  mode(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).PrevPd);
                        
                        
                        % Save MTF data
                        MTFdata = [MTFdata; PdData(idx,1).stimID  iseq  find(AMrates==this_rate) PdData(idx,2).(PlotThis) PdData(idx,1).(PlotThis) ];

%                         if strcmp(session,'MA') && channel==13
%                             keyboard
%                         end
                        
                        %% Add datapoints to population plots
                        
                        % CONTEXT COMPARISON
                        
                        % Each AM rate, Early/Late
                        figure(hfe); hold on
                        subplot(hsp(AMrates==this_rate,iseq)); hold on
                        
                        switch PlotThis
                            
                            case 'FR'                                
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,2).SEM,'horizontal','Color',colors(subMPH.ThisStimID(1),:))
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,1).SEM,'vertical','Color',colors(subMPH.ThisStimID(1),:))
                                scatter(PdData(idx,2).FR, PdData(idx,1).FR,dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                            case {'FF' 'TrV' 'VS' 'MP' 'CV' 'r_all' 'r_pdc'}
                                scatter(PdData(idx,2).(PlotThis), PdData(idx,1).(PlotThis),dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                        end
                        hold off
                        
                        % Each AM rate merged
                        figure(hfr); hold on
                        subplot(hspr(AMrates==this_rate)); hold on
                        
                        switch PlotThis
                            
                            case 'FR'
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,2).SEM,'horizontal','Color',colors(subMPH.ThisStimID(1),:))
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,1).SEM,'vertical','Color',colors(subMPH.ThisStimID(1),:))
                                scatter(PdData(idx,2).FR, PdData(idx,1).FR,dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                            case {'FF' 'TrV' 'VS' 'MP' 'CV' 'r_all' 'r_pdc'}
                                scatter(PdData(idx,2).(PlotThis), PdData(idx,1).(PlotThis),dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                        end
                        hold off
                        
                        % Overlay -- All
                        figure(hfo); hold on                        
                        switch PlotThis
                            
                            case 'FR'
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,2).SEM,'horizontal','Color',colors(subMPH.ThisStimID(1),:))
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,1).SEM,'vertical','Color',colors(subMPH.ThisStimID(1),:))
                                scatter(PdData(idx,2).FR, PdData(idx,1).FR,dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                            case {'FF' 'TrV' 'VS' 'MP' 'CV' 'r_all' 'r_pdc'}
                                scatter(PdData(idx,2).(PlotThis), PdData(idx,1).(PlotThis),dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)                                
                        end
                        hold off
                        
                        % Overlay -- Early/Late
                        figure(hfoel); hold on
                        subplot(hsel(iseq)); hold on
                        switch PlotThis
                            
                            case 'FR'
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,2).SEM,'horizontal','Color',colors(subMPH.ThisStimID(1),:))
                                errorbar(PdData(idx,2).FR, PdData(idx,1).FR, PdData(idx,1).SEM,'vertical','Color',colors(subMPH.ThisStimID(1),:))
                                scatter(PdData(idx,2).FR, PdData(idx,1).FR,dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                                
                            case {'FF' 'TrV' 'VS' 'MP' 'CV' 'r_all' 'r_pdc'}
                                scatter(PdData(idx,2).(PlotThis), PdData(idx,1).(PlotThis),dotsize,'o',...
                                    'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:),...
                                    'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
                        end
                        hold off
                        
                        
                        
                        % CONTEXT DIFF AS A FUNCTION OF STIMULUS
                        
                        if strcmp(select_stim,'all')
                            
                            diffData_unit = [diffData_unit; find(AMrates==this_rate) PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis) idx];
                            
                            % By AM rate
%                             figure(hfdr); hold on
%                             plot(find(AMrates==this_rate), PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis),...
%                                 'ok','MarkerSize',10)
%                             hold off
                            
                            hfdrData = [hfdrData; find(AMrates==this_rate) PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis)];
                            
                            % By distance to BMF-fr
                            if isfield(ThisResp,'iBMF_FR') && ~isempty(ThisResp.iBMF_FR)
                                
%                                 figure(hfdd); hold on
%                                 plot( find(AMrates==this_rate)-ThisResp.iBMF_FR, PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis),...
%                                     'ok','MarkerSize',10)
%                                 hold off
                                
                                hfddData = [hfddData; find(AMrates==this_rate)-ThisResp.iBMF_FR PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis)];
                            end
                        end
                        
                        
                        
                        % CONTEXT DIFF AS A FUNCTION OF PREVIOUS PD AND FR
                        
                        diffRaw  = PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis);
                        diffNorm = (PdData(idx,1).(PlotThis)-PdData(idx,2).(PlotThis))/PdData(idx,2).(PlotThis);
                        hfppData = [hfppData; this_rate  diffNorm  PdData(idx,1).prFR500  PdData(idx,1).prFR100  PdData(idx,1).prPd ];
                        
                        
                        
                        
                        %% Add datapoint to vectors for stats
                        
                        FinalDPs(this_rate==AMrates,iseq).Pdc = [ FinalDPs(this_rate==AMrates,iseq).Pdc   PdData(idx,2).(PlotThis)  ];
                        FinalDPs(this_rate==AMrates,iseq).IR  = [ FinalDPs(this_rate==AMrates,iseq).IR    PdData(idx,1).(PlotThis) ];
                        
                        DP_TrV(this_rate==AMrates,iseq).Pdc = [ DP_TrV(this_rate==AMrates,iseq).Pdc   PdData(idx,2).TrV ];
                        DP_TrV(this_rate==AMrates,iseq).IR  = [ DP_TrV(this_rate==AMrates,iseq).IR    PdData(idx,1).TrV ];
                        DP_FR(this_rate==AMrates,iseq).Pdc  = [ DP_FR(this_rate==AMrates,iseq).Pdc    PdData(idx,2).FR ];
                        DP_FR(this_rate==AMrates,iseq).IR   = [ DP_FR(this_rate==AMrates,iseq).IR     PdData(idx,1).FR ];
                        
                        
                        
                    end %iIRpd
                                        
                end  %this_rate
                
                
                %% ======================   MTF   ======================  
%                 figure(hfmtf); cla;  hold on
%                 for iirst = unique(MTFdata(:,1))'
%                     plot(MTFdata(MTFdata(:,2)==1 & MTFdata(:,1)==iirst,3), MTFdata(MTFdata(:,2)==1 & MTFdata(:,1)==iirst,4),'Color',[0 0 0],'LineWidth',2)
%                     plot(MTFdata(MTFdata(:,2)==2 & MTFdata(:,1)==iirst,3), MTFdata(MTFdata(:,2)==2 & MTFdata(:,1)==iirst,4),'Color',0.5*[1 1 1],'LineWidth',2)
%                     plot(MTFdata(MTFdata(:,2)==1 & MTFdata(:,1)==iirst,3), MTFdata(MTFdata(:,2)==1 & MTFdata(:,1)==iirst,5),'b','LineWidth',2)
%                     plot(MTFdata(MTFdata(:,2)==2 & MTFdata(:,1)==iirst,3), MTFdata(MTFdata(:,2)==2 & MTFdata(:,1)==iirst,5),'c','LineWidth',2)
%                 end
%                 keyboard
                
                
%                 % If want to plot lines between datapoints of same unit: 
%                 if strcmp(select_stim,'all')
%                     
%                     for ii = 1:max(diffData_unit(:,3))
%                         
%                         data = diffData_unit(diffData_unit(:,3)==ii,:);
%                         
%                         figure(hfdr); hold on
%                         plot(data(:,1), data(:,2), 'o-k','MarkerSize',10)
%                         hold off
%                         
%                         if isfield(ThisResp,'iBMF_FR') && ~isempty(ThisResp.iBMF_FR)
%                             
%                             figure(hfdd); hold on
%                             plot( data(:,1)-ThisResp.iBMF_FR, data(:,2), 'o-k','MarkerSize',10)
%                             hold off
%                         end 
%                         
%                     end
%                 end
                
                
                
                %%   End of datapoint
                                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects



%%  STATS

%% Individual plots -- AM rate and timing

figure(hfe); hold on
suptitle(sprintf('%s for each period\nselect stim: %s',PlotThis,select_stim))
try
for ir=1:5
    for is=1:2
        
        if numel(FinalDPs(ir,is).Pdc)>1
            
            % Effect of context
            
            % Get correlation of datapoints between contexts
            [r,p_r] = corrcoef( FinalDPs(ir,is).Pdc, FinalDPs(ir,is).IR);
            % Perform t-test of datapoints between contexts
            [~,p_t] = ttest2(   FinalDPs(ir,is).Pdc, FinalDPs(ir,is).IR);
            % Perform non-parametric (paired) Wilcoxon sign rank test of datapoints between contexts
            p_wsr   = signrank( FinalDPs(ir,is).Pdc, FinalDPs(ir,is).IR, 'method', 'exact');
            %         p_st =signtest(FinalDPs(ir,is).Pdc, FinalDPs(ir,is).IR); %similar
            %         to rank sign test but less powerful
            
            
            % Save stats
            FinalDPs(ir,is).stats.r     = r(1,2);
            FinalDPs(ir,is).stats.p_r   = p_r(1,2);
            FinalDPs(ir,is).stats.p_t   = p_t;
            FinalDPs(ir,is).stats.p_wsr = p_wsr;
            
        else
            % Save stats
            FinalDPs(ir,is).stats.r     = nan;
            FinalDPs(ir,is).stats.p_r   = nan;
            FinalDPs(ir,is).stats.p_t   = nan;
            FinalDPs(ir,is).stats.p_wsr = nan;
        end
        
        % Print in plots 
        subplot(hsp(ir,is)); hold on
%         title(sprintf('r=%2.2f, p=%2.2f\nwx sign rank p=%2.2f',r(1,2),p_r(1,2),p_wsr))
        title(sprintf('WSR p=%0.2e',p_wsr))
        
    end
    
    
    % Effect of timing
    
    if numel(FinalDPs(ir,1).Pdc)>1
        
        % Get correlation of datapoints from early and late times, of each context
        [r_Pdc,p_r_Pdc] =  corrcoef( FinalDPs(ir,1).Pdc, FinalDPs(ir,2).Pdc);
        [r_IR,p_r_IR]   =  corrcoef( FinalDPs(ir,1).IR,  FinalDPs(ir,2).IR);
        % Perform non-parametric (paired) Wilcoxon sign rank test of periodic/IR correlation coefficients
        p_wsr_Pdc       =  signrank( FinalDPs(ir,1).Pdc, FinalDPs(ir,2).Pdc, 'method', 'exact');
        p_wsr_IR        =  signrank( FinalDPs(ir,1).IR,  FinalDPs(ir,2).IR,  'method', 'exact');
        
        % Save stats
        FinalDPs(ir,1).stats.r_time_Pdc     =  r_Pdc(1,2);
        FinalDPs(ir,1).stats.r_time_IR      =  r_IR(1,2);
        FinalDPs(ir,1).stats.p_r_time_Pdc   =  p_r_Pdc(1,2);
        FinalDPs(ir,1).stats.p_r_time_IR    =  p_r_IR(1,2);
        FinalDPs(ir,1).stats.p_wsr_time_Pdc =  p_wsr_Pdc;
        FinalDPs(ir,1).stats.p_wsr_time_IR  =  p_wsr_IR;
        
    else
        % Save stats
        FinalDPs(ir,1).stats.r_time_Pdc     =  nan;
        FinalDPs(ir,1).stats.r_time_IR      =  nan;
        FinalDPs(ir,1).stats.p_r_time_Pdc   =  nan;
        FinalDPs(ir,1).stats.p_r_time_IR    =  nan;
        FinalDPs(ir,1).stats.p_wsr_time_Pdc =  nan;
        FinalDPs(ir,1).stats.p_wsr_time_IR  =  nan;
    end
    

    
end
hold off
catch
    keyboard
end



%% Overlay data - AM rates separate

figure(hfr); hold on
suptitle(sprintf('%s for each period\nselect stim: %s',PlotThis,select_stim))
try
    for ir=1:5
        
        if numel(FinalDPs(ir).Pdc)>1
            
            % Effect of context
            
            % Get correlation of datapoints between contexts
            [r,p_r] = corrcoef( FinalDPs(ir).Pdc, FinalDPs(ir).IR);
            % Perform t-test of datapoints between contexts
            [~,p_t] = ttest2(   FinalDPs(ir).Pdc, FinalDPs(ir).IR);
            % Perform non-parametric (paired) Wilcoxon sign rank test of datapoints between contexts
            p_wsr   = signrank( FinalDPs(ir).Pdc, FinalDPs(ir).IR, 'method', 'exact');
            
            % Save stats
            FinalDPs(ir,3).stats.r     = r(1,2);
            FinalDPs(ir,3).stats.p_r   = p_r(1,2);
            FinalDPs(ir,3).stats.p_t   = p_t;
            FinalDPs(ir,3).stats.p_wsr = p_wsr;
            
        else
            % Save stats
            FinalDPs(ir,3).stats.r     = nan;
            FinalDPs(ir,3).stats.p_r   = nan;
            FinalDPs(ir,3).stats.p_t   = nan;
            FinalDPs(ir,3).stats.p_wsr = nan;
        end
        
        % Print in plots
        subplot(hspr(ir)); hold on
        title(sprintf('WSR p=%0.2e, r=%0.2e p=%0.2e',p_wsr,r(1,2),p_r(1,2)))        
        
    end
    hold off
catch
    keyboard
end


%% Overlay data - Early/Late
p_wsr_Early = signrank([FinalDPs(:,1).Pdc],[FinalDPs(:,1).IR]);
p_wsr_Late  = signrank([FinalDPs(:,2).Pdc],[FinalDPs(:,2).IR]);
FinalDPs(1,1).stats.p_wsr_Early = p_wsr_Early;
FinalDPs(1,2).stats.p_wsr_Late  = p_wsr_Late;

figure(hfoel); hold on
subplot(hsel(1)); hold on
title(sprintf('%s for each %s period\nselect stim: %s\nWSR p=%0.2e',PlotThis,SeqPosStr{1},select_stim,p_wsr_Early))
subplot(hsel(2)); hold on
title(sprintf('%s for each %s period\nselect stim: %s\nWSR p=%0.2e',PlotThis,SeqPosStr{2},select_stim,p_wsr_Late))
hold off


%% Overlay data - All
p_wsr_ALL = signrank([FinalDPs.Pdc],[FinalDPs.IR]);
FinalDPs(1,1).stats.p_wsr_allRates = p_wsr_ALL;

figure(hfo);
title(sprintf('%s for each period\nselect stim: %s\nWSR p=%0.2e',PlotThis,select_stim,p_wsr_ALL))



%% Check relationship between TrV and FR
% 
% % Context diff then TrV vs FR
% figure;
% subplot(1,2,1); hold on
% plot(([DP_TrV.IR]-[DP_TrV.Pdc]),([DP_FR.IR]-[DP_FR.Pdc]),'xk' )
% xlabel('TrV(IR) - TrV(Pdc)')
% ylabel('FR(IR) - FR(Pdc)')
% [r,p]=corrcoef(([DP_TrV.IR]-[DP_TrV.Pdc]),([DP_FR.IR]-[DP_FR.Pdc]));
% title(sprintf('TrV context diff vs. FR context diff\nr = %0.3f, p = %0.2e',r(1,2),p(1,2)))
% 
% 
% % TrV/FR then Pdc vs IR
% normVar_Pdc = [DP_TrV.Pdc]./[DP_FR.Pdc];
% normVar_IR  = [DP_TrV.IR]./[DP_FR.IR];
% rm_i = find(isnan(normVar_Pdc));
% rm_i = [rm_i find(isnan(normVar_IR))];
% normVar_Pdc(rm_i) = [];
% normVar_IR(rm_i)  = [];
% 
% subplot(1,2,2); hold on
% plot(normVar_Pdc,normVar_IR,'xk' )
% xlabel('TrV/FR (Pdc)')
% ylabel('TrV/FR (IR)')
% [r,p]=corrcoef(normVar_Pdc,normVar_IR);
% title(sprintf('Normalized variability (TrV/FR), Pdc vs IR\nr = %0.3f, p = %0.2e',r(1,2),p(1,2)))



%% IR-Pdc diff against AM rate/distance to BMF
if strcmp(select_stim,'all')
    
    figure(hfdr); hold on
    plotSpread(hfdrData(:,2),'distributionIdx',hfdrData(:,1),'distributionColors','k')
    boxplot(hfdrData(:,2)',hfdrData(:,1)','Positions',unique(hfdrData(:,1))',...
        'Colors','r','notch','on','Symbol','','Labels',Info.stim_ID_key(2:6)')
    
%     distributionPlot(hfdrData(:,2),'group',hfdrData(:,1),'xNames',Info.stim_ID_key(2:6)','showMM',6)   
%     boxplot(hfdrData(:,2)',hfdrData(:,1)','Positions',unique(hfdrData(:,1))','Colors','r','notch','on')
%     h=findobj(gca,'Tag','Median'); j=findobj(gca,'Tag','Box');
%     for ih=1:numel(h)
%         h(ih).LineWidth = 4;
%         j(ih).LineWidth = 2;
%     end
    [p,ant,stats] = kruskalwallis(hfdrData(:,2)',hfdrData(:,1)','off');
    title(sprintf('%s context diff as a fct of AM rate\nKW p=%0.2e',PlotThis,p))
    hold off
    
    figure(hfdd); hold on
    
    plotSpread(hfddData(:,2),'distributionIdx',hfddData(:,1),'distributionColors','k')
    boxplot(hfddData(:,2)',hfddData(:,1)','Positions',unique(hfddData(:,1))','Colors','r','notch','on','Symbol','')

%     boxplot(hfddData(:,2)',hfddData(:,1)','Positions',unique(hfddData(:,1))','Colors','r','notch','on')
%     h=findobj(gca,'Tag','Median'); j=findobj(gca,'Tag','Box');
%     for ih=1:numel(h)
%         h(ih).LineWidth = 4;
%         j(ih).LineWidth = 2;
%     end
    [p,ant,stats] = kruskalwallis(hfddData(:,2)',hfddData(:,1)','off');
    title(sprintf('%s context diff as a fct of distance from BMF-fr\nKW p=%0.2e',PlotThis,p))
    hold off
end



%% Check for a relationship between IR/Pdc response and prev pd and prev FR
% 
% figure;
% hold on
% for ir=1:5 
%     subplot(1,5,ir); hold on
%     ihf = hfppData(:,1)==AMrates(ir);
%     plot(hfppData(ihf,3), hfppData(ihf,2), 'xk','LineWidth',2,'MarkerSize',10)
%     title([num2str(AMrates(ir)) ' Hz period'])
%     set(gca,'xlim',[0 35],'xscale','linear','ylim',2*[-1 1])
%     hold off
% end
% suptitle('prev 500 ms FR')
% hold off
% 
% 
% figure;
% hold on
% for ir=1:5 
%     subplot(1,5,ir); hold on
%     ihf = hfppData(:,1)==AMrates(ir);
%     plot(hfppData(ihf,4), hfppData(ihf,2), 'xk','LineWidth',2,'MarkerSize',10)
%     title([num2str(AMrates(ir)) ' Hz period'])
%     set(gca,'xlim',[0 35],'xscale','linear','ylim',2*[-1 1])
%     hold off
% end
% suptitle('prev 100 ms FR')
% hold off
% 
% figure;
% hold on
% for ir=1:5 
%     subplot(1,5,ir); hold on
%     ihf = hfppData(:,1)==AMrates(ir);
%     plot(hfppData(ihf,5), hfppData(ihf,2), 'xk','LineWidth',2,'MarkerSize',10)
%     title([num2str(AMrates(ir)) ' Hz period'])
%     set(gca,'xlim',[1 64],'xscale','log','ylim',2*[-1 1])
%     hold off
% end
% suptitle('prev Pd')
% hold off




%% Save figure(s)

savedir = fullfile(fn.processed,'MPHmetrics','lowFR');
if SUonly==1
    savedir = [savedir '/SU'];
end
if ~exist(savedir,'dir')
    mkdir([savedir '/svg']);
    mkdir([savedir '/eps']);
end
savename = sprintf('MPHs_Context_%s_%s',PlotThis,select_stim);

% eps
print_eps_kp(hfe,fullfile([savedir '/eps'],[savename '_each']))
print_eps_kp(hfr,fullfile([savedir '/eps'],[savename '_rates']))
print_eps_kp(hfo,fullfile([savedir '/eps'],[savename '_overlay']))
% svg
print_svg_kp(hfe,fullfile([savedir '/svg'],[savename '_each']))
print_svg_kp(hfr,fullfile([savedir '/svg'],[savename '_rates']))
print_svg_kp(hfo,fullfile([savedir '/svg'],[savename '_overlay']))

if strcmp(select_stim,'all')
    % eps
    print_eps_kp(hfdr,fullfile([savedir '/eps'],[savename '_indAMrate']))
    print_eps_kp(hfdd,fullfile([savedir '/eps'],[savename '_indBMFdiff']))
    % svg
    print_svg_kp(hfdr,fullfile([savedir '/svg'],[savename '_indAMrate']))
    print_svg_kp(hfdd,fullfile([savedir '/svg'],[savename '_indBMFdiff']))
end


% Also save stats
save(fullfile(savedir,savename),'FinalDPs','-v7.3')




end %function




