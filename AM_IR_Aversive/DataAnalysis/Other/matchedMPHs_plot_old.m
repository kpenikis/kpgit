function plotMPHs()
%
%  MPHanalyses( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%
%  KP, 2018-04, 2019-03
%



global fn
close  all

if nargin<1
    PlotThis = 'FR'; 'FF'; 'TrV'; 'VS'; 'RS';
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
elseif nargin==1
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
end

%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!



%% Set up figure

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
narrow = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
smallrect = [1 scrsz(4)/3 scrsz(3)/5 scrsz(4)/3];
halfscreen = [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2];


AMrates = [2 4 8 16 32];
ContextStr = {'Periodic' 'IR (AC)' 'IR (DB)'};
SeqPosStr = {'Early' 'Late'};

yminval = 0;
ymaxval = 70;

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

fn = set_paths_directories([],[],1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


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
        FinalDPs(ir,is).IR = [];
    end
end


%%
for iUn = 1:size(UnitInfo,1)
    
    %%% still must edit for merged units
    %     if numel(UnitInfo(iUn,:).Session{:})>2  %strncmp(UnitInfo.RespType{iUn},'merged',6)
    %         continue
    %     end
    
    subject = UnitInfo(iUn,:).Subject{:};
    session = UnitInfo(iUn,:).Session{:};
    channel = UnitData(iUn).Channel(1);
    clu     = UnitData(iUn).Clu(1);
    
    % Load data files
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % Get spiketimes and shift based on calculated integration time (KS)
    iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
    spiketimes = round(Clusters(iClu).spikeTimes*1000)';
    
    subjcol = [1 1 1];
        
    
    %%
    %~~~~~~~~~~~~~~~~~~~~~~~~
    % Collect FR over session
    %~~~~~~~~~~~~~~~~~~~~~~~~
    
    bs_hist = 1;
    bs_smth = 20;
    
    [Stream_FRsmooth,Stream_zscore,Stream_Spikes] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');
    
    
    % Check if unit has very low overall firing rate
    if mean(Stream_FRsmooth(TrialData.onset(1):end)) < 1
        continue
    end
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    N = N+1;
    
    
    %%
    
    % Get sound parameters
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    
    % Get all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    
    all_TDidx = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    allStim = unique(TrialData.trID(all_TDidx));
    
    
    
    %==========================================================
    % Set up empty data structure for this unit
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
                    
                    MPH_addrow = {istim pstid allPds(ipd) seqPos starttime Prev500msFR  {MPH_temp(ipd).raster} dBSPL LP};
                    MPH = [MPH; MPH_addrow];
                    
                end
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                clear MPH_temp
                
                
            end %pstid (prev stim id)
            
            
        end %ITIflag (Trial or ITI)
    end %istim
    
    MPH(1,:) = [];
    
    %                 mnNtrials = min(cell2mat(cellfun(@(x) size(x,1),MPH.raster,'UniformOutput',false)));
    
    
    
    %%  3) PLOT DATA
    
    
    for this_rate = AMrates
        
        % Set up
        Period = 1000/this_rate;
        theseIRs = unique(MPH.ThisStimID(MPH.ThisStimID>6))';
        PdData = struct;
        
        %==========================================================
        % Set up figure
        %==========================================================
        
        ir = find(AMrates==this_rate);
        
        hfe(ir)   =   figure;
        set(hfe(ir),'Position',halfscreen)
        
        for iseq=1:2
            for ic = (theseIRs-6)
                
                hsp(iseq,ic)=subplot(2,2,iseq*2-2+ic); hold on
                axis square
                
                if iseq==2
                    xlabel('Time in period')
                end
                
                if ic==1
                    ylabel(sprintf('%s\nSpikes/sec',SeqPosStr{iseq}))
                else
                    ylabel('Spikes/sec')
                end
            end
        end
        linkaxes(hsp,'xy')
        xlim([0 ceil(1000/this_rate)])
        ylim([yminval ymaxval])
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Get IR context first, and collect starttimes to match
        
        
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
                
                raster=[];
                
                for ipst = 1:numel(subMPH(subMPH.SeqPos==iseq,:).raster)
                    
                    raster = [raster; subMPH(subMPH.SeqPos==iseq,:).raster{ipst}];
                    
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
                        spktimes = [spktimes find(raster(it,:))];
                    end
                    VS(ipst) = vectorstrength(spktimes,Period);
                    
                end %ipst
                
                
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
                
                
                %================================
                %           PLOT IR
                %================================
                
                histMP = [];
                histbin = floor(size(raster,2)/30);
                histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
                
                
                subplot(hsp(iseq,IRstim-6)); hold on
                %                             title(sprintf('%s\n%i ms, %i tr',ContextStr{IRstim-5},PdData(idx,1).starttime, size(raster,1) ))
                titstrIR{iseq,IRstim-6} = sprintf('%s: %i ms, %i tr',ContextStr{IRstim-5},PdData(idx,1).starttime, size(raster,1) );
                
                bar(linspace(1,size(raster,2),length(histMP)),histMP,...
                    'FaceColor',colors(ir+1,:),'EdgeColor','none','BarWidth',1)
                
                hold off
                
                
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
            
            % Find matching starttime
            strtdiffs = PdData(idx,1).starttime - pdcstarttimes;
            [~,ipd] = min(strtdiffs(strtdiffs>=0));
            %                         [~,ipd] = min(abs(PdData(idx,1).starttime - pdcstarttimes));
            
            % Preallocate
            FR  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            SEM = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            FF  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            TrV = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            VS  = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            
            raster=[];
            
            for ipst = 1:numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster)
                
                raster = [raster; subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster{ipst}];
                
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
                    spktimes = [spktimes find(raster(it,:))];% + (it-1)*Period];
                end
                VS(ipst) = vectorstrength(spktimes,Period);
                
            end %ipst
            
            
            PdData(idx,2).stimID = this_rate;
            PdData(idx,2).seqPos = iseq;
            PdData(idx,2).starttime = mode(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Starttime);
            PdData(idx,2).FR  = mean(FR,'omitnan');
            PdData(idx,2).SEM = mean(SEM,'omitnan');
            PdData(idx,2).FF  = mean(FF,'omitnan');
            PdData(idx,2).TrV = mean(TrV,'omitnan');
            PdData(idx,2).VS  = mean(VS,'omitnan');
            
            
            
            %================================
            %           PLOT Pdc
            %================================
            
            histMP = [];
            histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
            
            %                         subplot(hsp(iseq,1)); hold on
            subplot(hsp(iseq,PdData(iIRpd,1).stimID-6));
            hold on
            
            %                         bar(linspace(1,size(raster,2),length(histMP)),histMP,...
            %                             'FaceColor','none','EdgeColor',0*[1 1 1],'BarWidth',1,'LineWidth',1)
            
            xx = linspace(1,size(raster,2),length(histMP));
            xx = reshape(mode(diff(xx))/2+[xx;xx],1,2*length(xx));
            if ~isempty(histMP)
                yy = reshape([histMP;histMP],1,2*length(histMP));
                yy(1) = []; yy(1) = 0; yy(end+1) = 0;
            else
                yy = zeros(size(xx));
            end
            
            plot(xx,yy,'-','Color',0.7*[1 1 1],'LineWidth',1.5)
            
            title(sprintf('%s\n%s: %i ms, %i tr',titstrIR{iseq,PdData(iIRpd,1).stimID-6},ContextStr{1},PdData(idx,2).starttime, size(raster,1) ))
            
            hold off
            
        end %iIRpd
        
        
        %% Finish and save plot
        
%         if isfield(ThisResp,'iBMF_FR') && isfield(ThisResp,'iBMF_VS')
%             suptitle( sprintf('%i Hz   |   %s %s %i %i   |   BMF-fr: %i Hz, BMF-vs: %i Hz\n ',this_rate,subject,session,channel,clu , AMrates(ThisResp.iBMF_FR), AMrates(ThisResp.iBMF_VS) ) )
%         elseif   isfield(ThisResp,'iBMF_FR')  &&  ~isfield(ThisResp,'iBMF_VS')
%             suptitle( sprintf('%i Hz   |   %s %s %i %i   |   BMF-fr: %i Hz, ns Sync\n ',this_rate,subject,session,channel,clu , AMrates(ThisResp.iBMF_FR) ) )
%         elseif  ~isfield(ThisResp,'iBMF_FR')  &&   isfield(ThisResp,'iBMF_VS')
%             suptitle( sprintf('%i Hz   |   %s %s %i %i   |   ns FR resp, BMF-vs: %i Hz\n ',this_rate,subject,session,channel,clu , AMrates(ThisResp.iBMF_VS) ) )
%         elseif  ~isfield(ThisResp,'iBMF_FR')  &&  ~isfield(ThisResp,'iBMF_VS')
%             suptitle( sprintf('%i Hz   |   %s %s %i %i   |   ns FR resp, ns Sync\n ',this_rate,subject,session,channel,clu ) )
%         end
%         
%         savedir = fullfile(fn.processed,'MPHplots','floor');
%         if SUonly==1
%             savedir = [savedir '/SU'];
%         end
%         if ~exist(savedir,'dir')
%             mkdir([savedir '/svg']);
%             mkdir([savedir '/eps']);
%         end
%         
%         savename = sprintf('MPH%iHz_%s_%s_ch%i_%i',this_rate,subject,session,channel,clu);
%         
%         print_eps_kp(hfe(ir),fullfile([savedir '/eps'],savename))
%         print_svg_kp(hfe(ir),fullfile([savedir '/svg'],savename))
        
        
    end  %this_rate
    close all
    
end %iUn



end %function




