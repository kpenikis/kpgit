function MPH_matchedComparison()
%
%  MPHanalyses( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%
%  KP, 2018-04, 2019-03
%



global fn AMrates rateVec_AC rateVec_DB trMin RateStream

%!!!!!!!!!!!!!!!!!
trMin   =  10;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!
AMrates = [2 4 8 16 32];
%!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

ContextStr = {'Periodic' 'IR (AC)' 'IR (DB)'};
SeqPosStr = {'Early' 'Late'};

USE_MEASURE = 'FR';
switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 0.01;
                axmax = 15*ceil(max([UnitData.BaseFR])/10);
            case 'z'
                axmin = -1;
                axmax = 1.5;
        end
    case {'TrV' 'FF'}
        units = ' ';
        axmin = 0.01;
        axmax = 4;
end



%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

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

% yvals = [5 10 20 30 50 75 100];
yvals = 2.^[0:7];


hf = figure;
set(hf,'Position',fullscreen)
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square


%% Preallocate

PairedData=[];

Pdata = table;
Pdata.Subject = ' ';
Pdata.Session = ' ';
Pdata.Channel = nan;
Pdata.Clu     = nan;
Pdata.IRstim  = nan;
Pdata.BaseFR  = nan;
Pdata.pval    = nan;
Pdata.predIR  = nan;
Pdata.obsIR   = nan;


for iUn = 1:numel(UnitData)
        
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    subjcol     = [1 1 1];
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
        
    
    % Load data files
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 || ~exist('TrialData','var')
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1 || ( ~exist('Spikes','var') || ~exist('Clusters','var') )
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if ~isfield(Info,'artifact')
        keyboard
    end
    
    % Get spiketimes and shift based on calculated integration time
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 + spkshift));  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 + spkshift)');
        
    end
    
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    % Set ymaxval
    yin = find( ( max(yvals, 2+3*max(nanmean(UnitData(iUn).FR_raw_tr,1)) ) - yvals )==0 );
    
    
    %% Get MPH data
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes,RateStream);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    %%  3) PLOT DATA
    
    for this_rate = AMrates(1:4)
        
        % Set up
        Period = 1000/this_rate;
        theseIRs = unique(MPH.ThisStimID(MPH.ThisStimID>6))';
        PdData = struct;
        
        
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
                
%                 histMP = [];
%                 histbin = floor(size(raster,2)/30);
%                 histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
%                 
%                 
%                 subplot(hsp(iseq,IRstim-6)); hold on
%                 %                             title(sprintf('%s\n%i ms, %i tr',ContextStr{IRstim-5},PdData(idx,1).starttime, size(raster,1) ))
%                 titstrIR{iseq,IRstim-6} = sprintf('%s: %i ms, %i tr',ContextStr{IRstim-5},PdData(idx,1).starttime, size(raster,1) );
%                 
%                 bar(linspace(1,size(raster,2),length(histMP)),histMP,...
%                     'FaceColor',colors(ir+1,:),'EdgeColor','none','BarWidth',1)
%                 
%                 hold off
%                 
                
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
            
%             histMP = [];
%             histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
%             
%             %                         subplot(hsp(iseq,1)); hold on
%             subplot(hsp(iseq,PdData(iIRpd,1).stimID-6));
%             hold on
%             
%             %                         bar(linspace(1,size(raster,2),length(histMP)),histMP,...
%             %                             'FaceColor','none','EdgeColor',0*[1 1 1],'BarWidth',1,'LineWidth',1)
%             
%             xx = linspace(1,size(raster,2),length(histMP));
%             xx = reshape(mode(diff(xx))/2+[xx;xx],1,2*length(xx));
%             if ~isempty(histMP)
%                 yy = reshape([histMP;histMP],1,2*length(histMP));
%                 yy(1) = []; yy(1) = 0; yy(end+1) = 0;
%             else
%                 yy = zeros(size(xx));
%             end
%             
%             plot(xx,yy,'-','Color',0.7*[1 1 1],'LineWidth',1.5)
%             
%             title(sprintf('%s\n%s: %i ms, %i tr',titstrIR{iseq,PdData(iIRpd,1).stimID-6},ContextStr{1},PdData(idx,2).starttime, size(raster,1) ))
%             
%             hold off
            
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
        
        
        %% Compare matched periods
        
        thesePds = find([PdData(:,1).seqPos]==2);
        
        for idx = thesePds
            
            if ~any(isnan([PdData(idx,:).FR]))
                
                plot(PdData(idx,2).FR, PdData(idx,1).FR, 'o','MarkerSize',15,'Color', colors(find(this_rate==AMrates)+1,:) )
            
                % Save values for stats later
                PairedData = [ PairedData; PdData(idx,2).FR PdData(idx,1).FR ];
            end
            
            % Should I compare against 'within periodic' again?
        end
        
        
        
    end  %this_rate
    
end %iUn


% Stats
pwsr    = signrank(PairedData(:,1),PairedData(:,2));
[rc,pc] = corrcoef(PairedData(:,1),PairedData(:,2));

title(sprintf('%s  |  wsr p=%0.2f, r=%0.2f p=%0.4e',USE_MEASURE,pwsr,rc(1,2),pc(1,2)))


% Save figure
savedir = fullfile(fn.figs,'MPHmatched');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,['matchedMPH_' USE_MEASURE]));



end %function




