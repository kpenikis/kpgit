function MPH_matchedComparison(USE_MEASURE)
%  MPH_matchedComparison
%
%   Gets all MPHs. Finds the closest matching Pdc period for each Irr
%   period. Also finds a different Pdc period to compare to this one. 
%   Option to plot MPH.
%
% KP, 2018-04, 2019-08
% 

% for adapting or facilitating cells, plot abs(Obs-Pred) as a function of 
% Pdc MP start time


% disp('Baseline-subtracted')

global fn AMrates trMin RateStream

%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!
colorswitch = 'PdcType'; 'AMrate'; 
%!!!!!!!!!!!!!!!!!
PlotMPH = 0;
%!!!!!!!!!!!!!!!!!

if nargin<1
    USE_MEASURE = 'FR';
end

%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

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


savedir = fullfile(fn.figs,'MPHmatched');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = ['matchedMPH_' USE_MEASURE ];


%% Figure settings

% set(0,'DefaultTextInterpreter','none')
% set(0,'DefaultAxesFontSize',14)
% 
% scrsz = get(0,'ScreenSize');  %[left bottom width height]
% fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
% 
% % Set colors
% colors = [ 250 250 250;...
%     84  24  69;...
%     120  10  41;...
%     181   0  52;...
%     255  87  51;...
%     255 153   0]./255;
% colors = [ colors; ...
%     [37  84 156]./255 ;...
%     [19 125 124]./255 ];
% 
% PlotMarkerSize = 15;
% 
% histbinedges = linspace(-sqrt(60^2 + 60^2)/2,sqrt(60^2 + 60^2)/2,40);


%% Preallocate

PD_PdcIR=[];
PD_IRIR=[];
PD_PdcPdc=[];

PdcRespType={};
PdcStartTime=[];
MPrate=[];
PRT = cell(numel(UnitData),5);
CellType={};


for iUn = 1:numel(UnitData)
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
    
    
    % Load data files
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 || ~exist('TrialData','var')
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info RateStream
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
%     end
%     if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1 || ( ~exist('Spikes','var') || ~exist('Clusters','var') )
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if ~isfield(Info,'artifact') || ~exist('RateStream','var')
        keyboard
    end
    
    % Get spiketimes and shift based on calculated integration time
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift));  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
        
    end
    
    
    if isempty(UnitData(iUn).iBMF_FR)
        continue
    end
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    
    
    %% Get MPH data
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes,RateStream);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    
    
    %%  3) GET MATCHED PERIODS, COLLECT DATA
    
    for this_rate = AMrates %(UnitData(iUn).iBMF_FR)
        
        % Set up
        Period    = 1000/this_rate;
        PdcIRData = struct;
        PdcData   = struct;
        IRData    = struct;
        T         = 250;
        alfa      = 0.05;
        
        % Get Pdc response type                      [ early  late  p_val ]
        PRT{iUn,this_rate==AMrates} = 'n';
        if ~isempty(UnitData(iUn).DeltaNspk) && ~isempty(UnitData(iUn).DeltaNspk{this_rate==AMrates})
            if (UnitData(iUn).DeltaNspk{this_rate==AMrates}(3))<alfa && (UnitData(iUn).DeltaNspk{this_rate==AMrates}(1) > UnitData(iUn).DeltaNspk{this_rate==AMrates}(2))
                % adapting
                PRT{iUn,this_rate==AMrates} = 'A';
            elseif (UnitData(iUn).DeltaNspk{this_rate==AMrates}(3))<alfa && (UnitData(iUn).DeltaNspk{this_rate==AMrates}(1) < UnitData(iUn).DeltaNspk{this_rate==AMrates}(2))
                % facilitating
                PRT{iUn,this_rate==AMrates} = 'F';
            else
                PRT{iUn,this_rate==AMrates} = 'S';
            end
        end
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Get instances of this period in IRREGULAR context first
        
        theseIRs = unique(MPH.ThisStimID(MPH.ThisStimID>6))';
        
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
                
                
                % Average over previous stimuli and save for later
                idx = 2*(IRstim-7) + iseq;
                PdcIRData(idx,1).stimID    = IRstim;
                PdcIRData(idx,1).seqPos    = iseq;
                PdcIRData(idx,1).starttime = mode(subMPH(subMPH.SeqPos==iseq,:).Starttime);
                PdcIRData(idx,1).FR        = mean(FR,'omitnan');% - UnitData(iUn).BaseFR;
                PdcIRData(idx,1).SEM       = mean(SEM,'omitnan');
                PdcIRData(idx,1).FF        = mean(FF,'omitnan');
                PdcIRData(idx,1).TrV       = mean(TrV,'omitnan');
                PdcIRData(idx,1).VS        = mean(VS,'omitnan');
                
                
                %% PLOT IR MPH
                %================================
                %           PLOT IR
                %================================
                
                if PlotMPH && idx==2
                    hfmph=figure; hold on
                    
                    histMP = [];
                    histbin = floor(size(raster,2)/30);
                    histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
                    
                    subplot(1,3,1); 
                    bar(linspace(1,size(raster,2),length(histMP)),histMP,...
                        'FaceColor','k','EdgeColor','none','BarWidth',1)
                    axis square
                    set(gca,'xtick',[],'ytick',[],'Color','none')
                    box on
                    xlim([0 size(raster,2)])
                    ylim([0 40])
                    title('IR')
                    xlabel(sprintf('%i ms',PdcIRData(idx,1).starttime))
                    
                end
                
            end %iseq
            
        end %IRstim -- IRREGULAR first
        
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Now get matching PERIODIC context data
        
        subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate,:);
        pdcstarttimes = unique(subMPH.Starttime)';
        
        % SKIP IF TOO CLOSE TO TRIAL TRANSITION
        pdcstarttimes(pdcstarttimes<T) = nan;
        
        for iIRpd = 1:size(PdcIRData,1)
            
            if isempty(PdcIRData(iIRpd,1).stimID)
                continue
            end
            
            IRstim = PdcIRData(iIRpd,1).stimID;
            iseq = PdcIRData(iIRpd,1).seqPos;
            idx = 2*(IRstim-7) + iseq;
            if idx~=iIRpd, keyboard, end
            
            % Find matching starttime
            strtdiffs = abs(PdcIRData(idx,1).starttime - pdcstarttimes);
            [STdiff,ipd] = min(strtdiffs);
            if STdiff>1000
                aaa=234;
            end
            
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
            
            PdcIRData(idx,2).stimID    = find(AMrates==this_rate)+1;
            PdcIRData(idx,2).seqPos    = iseq;
            PdcIRData(idx,2).starttime = mode(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Starttime);
            PdcIRData(idx,2).FR        = mean(FR,'omitnan');% - UnitData(iUn).BaseFR;
            PdcIRData(idx,2).SEM       = mean(SEM,'omitnan');
            PdcIRData(idx,2).FF        = mean(FF,'omitnan');
            PdcIRData(idx,2).TrV       = mean(TrV,'omitnan');
            PdcIRData(idx,2).VS        = mean(VS,'omitnan');
            
            
            %================================
            %           PLOT Pdc
            %================================
            if PlotMPH && idx==2
                
                histMP = [];
                histbin = floor(size(raster,2)/30);
                histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
                
                subplot(1,3,2);
                bar(linspace(1,size(raster,2),length(histMP)),histMP,...
                    'FaceColor','k','EdgeColor','none','BarWidth',1)
                axis square
                xlim([0 size(raster,2)])
                ylim([0 40])
                set(gca,'xtick',[],'ytick',[],'Color','none')
                box on
                title('Pdc 1')
                xlabel(sprintf('%i ms',PdcIRData(idx,2).starttime))
                
            end
            
            
            %---------------------------------------
            %  Also get MPH for Pdc-Pdc comparison
            %---------------------------------------
            
            % 2 Hz: always compare to first period, but only those preceded by 4 Hz
            if this_rate==2
                subsubMPH = subMPH(subMPH.Starttime==0,:);
                subsubMPH = subsubMPH(subsubMPH.PrevStimID==3,:);
            
            % 4 Hz: allow first period to be chosen, when preceded by 4 Hz
            elseif this_rate==4
                nt=0; int=0;
                while nt<trMin
                    int=int+1;
                    ipd2s = 1:numel(pdcstarttimes);
                    ipd2s(ipd2s==ipd) = [];
                    ipd2 = ipd2s(randperm(numel(ipd2s),1));
                    if ipd2==1
                        subsubMPH = subMPH(subMPH.Starttime==0,:);
                        subsubMPH = subsubMPH(subsubMPH.PrevStimID==3,:);
                    else
                        subsubMPH = subMPH(subMPH.Starttime==pdcstarttimes(ipd2),:);
                    end
                    nt = sum(cellfun(@(x) size(x,1),subsubMPH.raster));
                    if int>10
                        warning('issue finding an appropriate MPH to compare')
                        keyboard
                        break
                    end
                end
                
            % 8,16,32 Hz: only choose from periods after the first 250 ms
            else
                ipd2s = 1:numel(pdcstarttimes);
                ipd2s(ipd2s==ipd) = [];
                ipd2s(isnan(pdcstarttimes)) = [];
                ipd2 = ipd2s(randperm(numel(ipd2s),1));
                subsubMPH = subMPH(subMPH.Starttime==pdcstarttimes(ipd2),:);
            end
            
            if isempty(subsubMPH) || isempty(subMPH)
                keyboard
            end
            
            % Preallocate
            FR  = nan(1,numel(subsubMPH.raster));
            SEM = nan(1,numel(subsubMPH.raster));
            FF  = nan(1,numel(subsubMPH.raster));
            TrV = nan(1,numel(subsubMPH.raster));
            VS  = nan(1,numel(subsubMPH.raster));
            
            raster=[];
            
            for ipst = 1:numel(subsubMPH.raster)
                
                raster = [raster; subsubMPH.raster{ipst}];
                
                if ~(ceil(Period)==size(raster,2))
                    keyboard
                end
                
                FR(ipst)   = sum(sum( raster ))/size(raster,1)/Period*1000;
                SEM(ipst)  = std(sum( raster,2)/Period*1000) / sqrt(size(raster,1));
                FF(ipst)   = var(sum( raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                TrV(ipst)  = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                
                % Reformat raster for VS calculation
                spktimes=[];
                for it = 1:size(raster,1)
                    spktimes = [spktimes find(raster(it,:))];% + (it-1)*Period];
                end
                VS(ipst) = vectorstrength(spktimes,Period);
                
            end %ipst
            
            PdcData(idx,1).stimID    = find(AMrates==this_rate)+1;
            PdcData(idx,1).seqPos    = iseq;
            PdcData(idx,1).starttime = mode(subsubMPH.Starttime);
            PdcData(idx,1).FR        = mean(FR,'omitnan');% - UnitData(iUn).BaseFR;
            PdcData(idx,1).SEM       = mean(SEM,'omitnan');
            PdcData(idx,1).FF        = mean(FF,'omitnan');
            PdcData(idx,1).TrV       = mean(TrV,'omitnan');
            PdcData(idx,1).VS        = mean(VS,'omitnan');
            
            
            %================================
            %          PLOT Pdc 2
            %================================
            if PlotMPH && idx==2
                
                histMP = [];
                histbin = floor(size(raster,2)/30);
                histMP = binspikecounts(mean(raster,1),histbin)/histbin*1000;
                
                subplot(1,3,3);
                bar(linspace(1,size(raster,2),length(histMP)),histMP,...
                    'FaceColor','k','EdgeColor','none','BarWidth',1)
                axis square
                xlim([0 size(raster,2)])
                ylim([0 40])
                set(gca,'xtick',[],'ytick',[],'Color','none')
                box on
                title('Pdc 2')
                xlabel(sprintf('%i ms',PdcData(idx,1).starttime))
                
                keyboard
                % Finish and save MPH plot
                print_eps_kp(gcf,fullfile(savedir,'MPHplots_iUn84_4hz_7_2'));
                
            end
            
        end %iIRpd
        
        
        
        
        %  =========================================================
        %% ==============   Compare matched periods   ==============
        %  =========================================================
        
%         thesePds = find([PdcIRData(:,1).seqPos]==2);
        thesePds = unique([PdcIRData(:,1).seqPos]);
        
        for idx = 1:size(PdcIRData,1)
            
            if ~any(isnan([PdcIRData(idx,:).FR])) && ~isempty([PdcIRData(idx,:).FR])
                
                PdcRespType{end+1,1} =  PRT{iUn,this_rate==AMrates} ;
                PdcStartTime(end+1)  =  PdcIRData(idx,2).starttime ;
                MPrate(end+1)        =  find(this_rate==AMrates) ;
                if UnitInfo(iUn,:).TroughPeak<0.5
                    CellType{end+1,1} = 'N';
                else
                    CellType{end+1,1} = 'B';
                end

                
                %---------
                % Pdc-Pdc
%                 subplot(hs(2)); hold on
%                 plot(PdcIRData(idx,2).(USE_MEASURE), PdcData(idx,1).(USE_MEASURE), 'o','MarkerSize',15,'Color','k' )
                
                PD_PdcPdc  = [ PD_PdcPdc; PdcIRData(idx,2).(USE_MEASURE) PdcData(idx,1).(USE_MEASURE) ];
                
                %---------
                % Pdc-IR
%                 subplot(hs(1)); hold on
%                 plot(PdcIRData(idx,2).(USE_MEASURE), PdcIRData(idx,1).(USE_MEASURE), 'o','MarkerSize',15,'Color',plotcol )
                
                PD_PdcIR  = [ PD_PdcIR; PdcIRData(idx,2).(USE_MEASURE) PdcIRData(idx,1).(USE_MEASURE) ];
                
                %---------
                % IR-IR
                if idx<size(PdcIRData,1) && ~isempty(PdcIRData(idx+1,1).(USE_MEASURE))
%                     subplot(hs(3)); hold on
%                     plot(PdcIRData(idx+1,1).(USE_MEASURE), PdcIRData(idx,1).(USE_MEASURE), 'o','MarkerSize',15,'Color','k' )
                    
                    PD_IRIR  = [ PD_IRIR; PdcIRData(idx+1,1).(USE_MEASURE) PdcIRData(idx,1).(USE_MEASURE) ];
                elseif idx==size(PdcIRData,1) && size(PdcIRData,1)>2 && ~isempty(PdcIRData(1,1).(USE_MEASURE))
%                     subplot(hs(3)); hold on
%                     plot(PdcIRData(1,1).(USE_MEASURE), PdcIRData(idx,1).(USE_MEASURE), 'o','MarkerSize',15,'Color','k' )
                    
                    PD_IRIR  = [ PD_IRIR; PdcIRData(1,1).(USE_MEASURE) PdcIRData(idx,1).(USE_MEASURE) ];
                end %IR-IR
            end %if ~empty
        end  %each Irr MPH
        
        
    end  %this_rate
    
end %iUn

%% 

% Save data
save(fullfile(savedir,savename),'PD_PdcIR','PD_PdcPdc','PD_IRIR','PdcRespType','PdcStartTime','MPrate','CellType','-v7.3')

return

% Load already saved data
% load(fullfile(savedir,'matchedMPH_FR'))


%%

% Filter data by AM rate

% load(fullfile(savedir,'matchedMPH_FR'))
% 
% irate = 5;
% 
% PD_PdcIR     = PD_PdcIR(MPrate==irate,:);
% PD_PdcPdc    = PD_PdcPdc(MPrate==irate,:);
% % PD_IRIR      = PD_IRIR(MPrate==irate,:);
% PdcRespType  = PdcRespType(MPrate==irate);
% PdcStartTime = PdcStartTime(MPrate==irate);
% 
% figure; 
% plot(diff(PD_PdcPdc),diff(PD_PdcIR),'ko')


%% MAIN PLOTS

%::::::::::::::::::
%     SET UP
%::::::::::::::::::
hf = figure;
set(hf,'Position',fullscreen,'NextPlot','add')

hs(1)=subplot(6,4,[1 5 9]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-IR')

% Pdc-Pdc
hs(2)=subplot(6,4,[2 6 10]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-Pdc')

% Irr-Irr  
hs(3)=subplot(6,4,[3 7 11]);
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('IR-IR')

% Pdc-IR density plot
hs(4)=subplot(6,4,[13 17 21]);
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-IR')

% Pdc-Pdc density plot
hs(5)=subplot(6,4,[14 18 22]);
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('Pdc-Pdc')

% IR-IR density plot   OR   starttimeXdiff 
hs(6)=subplot(6,4,[15 19 23]);
set(gca,'Color','none','tickdir','out','ticklength',[0.025 0.025])
axis square
title('IR-IR')


% DIFFERENCE HISTOGRAMS

hs(7)=subplot(6,4,[4 8]);
hold on
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
set(gca,'Color','none','xtick',[],'ytick',[])
title('Pdc-IR')

hs(8)=subplot(6,4,[12 16]);
hold on
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
set(gca,'Color','none','xtick',[],'ytick',[])
title('Pdc-Pdc')

% diff histogram by 
hs(9)=subplot(6,4,[20 24]);
hold on
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
set(gca,'Color','none','xtick',[],'ytick',[])
title('IR-IR')

suptitle([num2str(AMrates(irate)) ' Hz'])


%% :::::::::::::::::
%     PLOT DATA  
%  :::::::::::::::::

%-----------------------
% Pdc-IR
subplot(hs(1)); hold on
% plot(PD_PdcIR(:,1), PD_PdcIR(:,2),'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor','k','MarkerEdgeColor','none')
plot(PD_PdcIR(strcmp(PdcRespType,'n'),1), PD_PdcIR(strcmp(PdcRespType,'n'),2), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.4*[1 1 0.9],'MarkerEdgeColor','none');
plot(PD_PdcIR(strcmp(PdcRespType,'S'),1), PD_PdcIR(strcmp(PdcRespType,'S'),2),  'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.01*[1 0.9 0.9],'MarkerEdgeColor','none');
plot(PD_PdcIR(strcmp(PdcRespType,'A'),1), PD_PdcIR(strcmp(PdcRespType,'A'),2),  'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[1 0.647 0],'MarkerEdgeColor','none');
plot(PD_PdcIR(strcmp(PdcRespType,'F'),1), PD_PdcIR(strcmp(PdcRespType,'F'),2),  'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[0 0.6 0],'MarkerEdgeColor','none');

%--Diff hist
subplot(hs(7)); hold on
histogram(PD_PdcIR(:,1)-PD_PdcIR(:,2),histbinedges,'Normalization','pdf','FaceColor','k','EdgeColor','none','FaceAlpha',1);
histogram(PD_PdcPdc(:,1)-PD_PdcPdc(:,2),histbinedges,'Normalization','pdf','FaceColor','none','EdgeColor','c','FaceAlpha',1,'DisplayStyle','stairs','LineWidth',2);
plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
ylim([0 0.3])

%--2d hist plot
subplot(hs(4)); hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
ihst(1)=histogram2(PD_PdcIR(:,1),PD_PdcIR(:,2),PlotMarkerSize,'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none');
cmocean('-gray')
% clims = get(gca,'CLim');
% set(gca,'CLim',clims+[-clims(2)*0.1 0])
set(gca,'CLim',[-10 100])


%-----------------------
% Pdc-Pdc
subplot(hs(2)); hold on
plot(PD_PdcPdc(:,1), PD_PdcPdc(:,2), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.01*[1 0.9 0.9],'MarkerEdgeColor','none');

%--Diff hist
subplot(hs(8)); hold on
histogram(PD_PdcPdc(:,1)-PD_PdcPdc(:,2),histbinedges,'Normalization','pdf','FaceColor','k','EdgeColor','none','FaceAlpha',1);
histogram(PD_PdcIR(:,1)-PD_PdcIR(:,2),histbinedges,'Normalization','pdf','FaceColor','none','EdgeColor','r','FaceAlpha',1,'DisplayStyle','stairs','LineWidth',2);
plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
ylim([0 0.3])

%--2d hist plot
% subplot(hs(5)); hold on
% plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
% ihst(2)=histogram2(PD_PdcPdc(:,1),PD_PdcPdc(:,2),PlotMarkerSize,'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none');
% cmocean('-gray')
% set(gca,'CLim',[-10 100])
% set(gca,'CLim',clims+[-clims(2)*0.1 0])


%-----------------------
% IR-IR
subplot(hs(3)); hold on
plot(PD_IRIR(:,1), PD_IRIR(:,2), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',0.01*[1 0.9 0.9],'MarkerEdgeColor','none');

% %--Diff hist
% subplot(hs(9)); hold on
% histogram(PD_IRIR(:,1)-PD_IRIR(:,2),histbinedges,'Normalization','pdf','FaceColor','k','EdgeColor','none','FaceAlpha',1);
% plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
% ylim([0 0.3])
% 
% %--2d hist plot
% subplot(hs(6)); hold on
% plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
% ihst(3)=histogram2(PD_IRIR(:,1),PD_IRIR(:,2),PlotMarkerSize,'DisplayStyle','tile','ShowEmptyBins','off','EdgeColor','none');
% cmocean('-gray')
% set(gca,'CLim',[-10 100])


%-------------------------------------------
% REPLACE density plots with:
% Difference as a function of start time 
%-------------------------------------------

Sdiffs = PD_PdcIR(strcmp(PdcRespType,'S'),1) - PD_PdcIR(strcmp(PdcRespType,'S'),2);
Fdiffs = PD_PdcIR(strcmp(PdcRespType,'F'),1) - PD_PdcIR(strcmp(PdcRespType,'F'),2);
Adiffs = PD_PdcIR(strcmp(PdcRespType,'A'),1) - PD_PdcIR(strcmp(PdcRespType,'A'),2);
ndiffs = PD_PdcIR(strcmp(PdcRespType,'n'),1) - PD_PdcIR(strcmp(PdcRespType,'n'),2);

subplot(hs(5)); hold on
plot([0 1000],[0 0],'Color',[0.7 0.7 0.7])
plot(PdcStartTime(strcmp(PdcRespType,'S')), Sdiffs, 'o','MarkerSize',PlotMarkerSize/2,'MarkerFaceColor',0.01*[1 0.9 0.9],'MarkerEdgeColor','none');
plot(PdcStartTime(strcmp(PdcRespType,'F')), Fdiffs, 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[0 0.6 0],'MarkerEdgeColor','none');
plot(PdcStartTime(strcmp(PdcRespType,'A')), Adiffs, 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[1 0.647 0],'MarkerEdgeColor','none');
xlim([0 1000])
ylim(40*[-1 1])
title('start time X diff')

subplot(hs(6)); hold on
plot([0 1000],[0 0],'Color',[0.7 0.7 0.7])
plot(PdcStartTime(strcmp(PdcRespType,'S')), Sdiffs / (PD_PdcIR(strcmp(PdcRespType,'S'),1) + PD_PdcIR(strcmp(PdcRespType,'S'),2)), 'o','MarkerSize',PlotMarkerSize/2,'MarkerFaceColor',0.01*[1 0.9 0.9],'MarkerEdgeColor','none');
plot(PdcStartTime(strcmp(PdcRespType,'F')), Fdiffs / (PD_PdcIR(strcmp(PdcRespType,'F'),1) + PD_PdcIR(strcmp(PdcRespType,'F'),2)), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[0 0.6 0],'MarkerEdgeColor','none');
plot(PdcStartTime(strcmp(PdcRespType,'A')), Adiffs / (PD_PdcIR(strcmp(PdcRespType,'A'),1) + PD_PdcIR(strcmp(PdcRespType,'A'),2)), 'o','MarkerSize',PlotMarkerSize,'MarkerFaceColor',[1 0.647 0],'MarkerEdgeColor','none');
xlim([0 1000])
ylim(0.8*[-1 1])
title('start time X norm diff')

subplot(hs(9)); hold on
histogram(Sdiffs, histbinedges, 'Normalization','pdf','FaceColor',0.01*[1 0.9 0.9],'EdgeColor','none','FaceAlpha',1); 
histogram(Fdiffs, histbinedges, 'Normalization','pdf','FaceColor', [0 0.6 0],'EdgeColor','none','FaceAlpha',1); 
histogram(Adiffs, histbinedges, 'Normalization','pdf','FaceColor', [1 0.647 0],'EdgeColor','none','FaceAlpha',1);
plot([0 0],[0 0.3],'Color',[0.7 0.7 0.7])
ylim([0 0.3])
xlim(sqrt(60^2 + 60^2)/2*[-1 1])
title('Diffs by Pdc resp type')



%% :::::::::::::::::
%     STATISTICS  
%  :::::::::::::::::

% inan = isnan(PD_PdcIR(:,1).*PD_PdcIR(:,2));
% PD_PdcIR(inan,:)  = [];
% PdcRespType(inan) = [];

% Stats                      [ Pdc  IR ]
[rc,pc] = corrcoef(PD_PdcIR(:,1), PD_PdcIR(:,2));
pwsr    = signrank(PD_PdcIR(:,1), PD_PdcIR(:,2));
pwrs    = ranksum(PD_PdcIR(:,1),  PD_PdcIR(:,2));
[~,ppt] = ttest(PD_PdcIR(:,1),    PD_PdcIR(:,2));
[~,pt2] = ttest2(PD_PdcIR(:,1),   PD_PdcIR(:,2));


subplot(hs(1)); hold on
stattext = sprintf('WSR p=%0.2f',pwsr);
text(axmax/4*3,axmax/15,stattext)
% if pc(1,2)<0.0000000001
%     title(sprintf('%s: Pdc vs IR   |   wsr p=%0.2f, r=%0.2f p<<0.0001',USE_MEASURE,pwsr,rc(1,2)))
% else
%     title(sprintf('%s: Pdc vs IR   |   wsr p=%0.2f, r=%0.2f p=%0.4e',USE_MEASURE,pwsr,rc(1,2),pc(1,2)))
% end


% Stats                      [ Pdc  Pdc ]
[rc,pc] = corrcoef(PD_PdcPdc(:,1), PD_PdcPdc(:,2));
pwsr    = signrank(PD_PdcPdc(:,1), PD_PdcPdc(:,2));

subplot(hs(2)); hold on
stattext = sprintf('WSR p=%0.2f',pwsr);
text(axmax/4*3,axmax/15,stattext)
% if pc(1,2)<0.0000000001
%     title(sprintf('%s: Pdc vs Pdc   |   wsr p=%0.2f, r=%0.2f p<<0.0001',USE_MEASURE,pwsr,rc(1,2)))
% else
%     title(sprintf('%s: Pdc vs Pdc   |   wsr p=%0.2f, r=%0.2f p=%0.4e',USE_MEASURE,pwsr,rc(1,2),pc(1,2)))
% end


% Stats                       [ IR  IR ]
[rc,pc] = corrcoef(PD_IRIR(:,1), PD_IRIR(:,2));
pwsr    = signrank(PD_IRIR(:,1), PD_IRIR(:,2));

subplot(hs(3)); hold on
stattext = sprintf('WSR p=%0.2f',pwsr);
text(axmax/4*3,axmax/15,stattext)
% if pc(1,2)<0.0000000001
%     title(sprintf('%s: IR vs IR   |   wsr p=%0.2f, r=%0.2f p<<0.0001',USE_MEASURE,pwsr,rc(1,2)))
% else
%     title(sprintf('%s: IR vs IR   |   wsr p=%0.2f, r=%0.2f p=%0.4e',USE_MEASURE,pwsr,rc(1,2),pc(1,2))) 
% end


%~~~~~~~~~~~~~~
% DIFFERENCES
%~~~~~~~~~~~~~~

% Pdc-Irr vs Pdc-Pdc

PdcDiffs  = PD_PdcPdc(:,1)-PD_PdcPdc(:,2);
PdIRDiffs = PD_PdcIR(:,1)-PD_PdcIR(:,2);

% Normality test
clear h_l
h_l(1) = lillietest(PdcDiffs);  % h=1 means distr NOT normal
h_l(2) = lillietest(PdIRDiffs);
if any(h_l==1)% non-parametric
    p_rs       = ranksum(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','BrownForsythe','display','off');
    stattext   = sprintf('ranksum p=%0.2f \nBrownForsythe p=%0.2e',p_rs,p_var);
else          % parametric
    [~,p_t]    = ttest2(PdcDiffs,PdIRDiffs);
    p_var      = vartestn([PdcDiffs PdIRDiffs],'TestType','Bartlett','display','off');
    [~,p_var]  = vartest2(PdcDiffs, PdIRDiffs);
    stattext   = sprintf('ttest p=%0.2f \nF test p=%0.2e',p_t,p_var);
    var_PdcPdc = var(PdcDiffs);
    var_PdcIR  = var(PdIRDiffs);
end

subplot(hs(8)); hold on
text(-40,0.28,stattext)


% Pdc-Irr by Pdc response type

[C,IA,IC]=unique(PdcRespType,'stable');

% Normality test
clear h_l
h_l(1) = lillietest(Sdiffs);   % h=1 means distr NOT normal
h_l(2) = lillietest(Fdiffs);
h_l(3) = lillietest(Adiffs);
h_l(4) = lillietest(ndiffs);

if any(h_l==1) % non-parametric
    [p_grp,tbl,stats] = kruskalwallis(PdIRDiffs,IC);
    grpstatstr = sprintf('kw test p=%0.2e',p_grp);
else           % parametric
    [p_grp,tbl,stats] = anovan(PdIRDiffs,IC);
    grpstatstr = sprintf('anova p=%0.2e',p_grp);
end
% multcompare(stats)
%     [h_tf,p_tf]=ttest2(Sdiffs(round(linspace(1,length(Sdiffs),numel(Fdiffs)))),Fdiffs);
%     [h_ta,p_ta]=ttest2(Sdiffs(round(linspace(1,length(Sdiffs),numel(Adiffs)))),Adiffs);

if h_l(3)==1
    avgA = median(Adiffs);
else
    avgA = mean(Adiffs);
end

figure(hf); hold on
subplot(hs(9)); hold on
text(-40,0.28,grpstatstr)
plot(avgA,0,'^k')


% Difference correlation with Pdc start time

[corr_A,pc_A] = corr(PdcStartTime(strcmp(PdcRespType,'A'))', Adiffs);

subplot(hs(5)); hold on
text(25,35,sprintf('A: r=%0.2f, p=%0.2f',corr_A,pc_A))


%% SAVE FIG AND DATA

% Save figure
print_eps_kp(hf,fullfile(savedir,['matchedMPH_' USE_MEASURE '_' num2str(AMrates(irate)) 'hz' ]));

keyboard
end %function




