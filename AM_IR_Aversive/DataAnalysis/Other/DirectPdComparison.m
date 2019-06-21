function DirectPdComparison( PlotThis, select_stim, UnitFilt )
%
%  analyzeMPH_withinContext( PlotThis, select_stim )
%   Control for analyzeMPH, which compares metrics of period responses
%   across Pdc/IR contexts. 
%
%  KP, 2018-05
%

global fn AMrates rateVec_AC rateVec_DB trMin
close  all
rng('shuffle')

if nargin<1
    PlotThis = 'covSpk'; 'FR'; 'FF'; 'TrV'; 'VS'; 'MP'; 'CrSc'; 'crW'; 'Rel'; 'CV'; 'OWR'; 'r_all'; 'r_pdc';
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
    UnitFilt =    'sparse';  'sustained'; 'gap'; 'merged'; 'na';  
elseif nargin==1
    select_stim = 'iBMF_FR'; 'iBMF_VS'; 'iSync'; 'all';
    UnitFilt =    'sparse';  'sustained'; 'gap'; 'merged'; 'na';  
elseif nargin==2
    UnitFilt =    'sparse';  'sustained'; 'gap'; 'merged'; 'na';  
end


%!!!!!!!!!!!!!!!!!
AMrespFilt =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
trMin   =  9;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!
OWbound = pi/3;
%!!!!!!!!!!!!!!!!!



%% Load Unit files

fn = set_paths_directories;

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
        FinalDPs(ir,is).IR  = [];
        DP_TrV(ir,is).Pdc = [];
        DP_TrV(ir,is).IR  = [];
        DP_FR(ir,is).Pdc  = [];
        DP_FR(ir,is).IR   = [];
    end
end


% Empty data matrices
Pdcmat  =[];
IRmat   =[];
PdIRmat =[];
hfdrData=[];
hfddData=[];
HistoryData = [];



%% Set options


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
narrow     = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
smallrect  = [1 scrsz(4)/3 scrsz(3)/5 scrsz(4)/3];
widerect = [1 scrsz(4)/3 scrsz(3) scrsz(4)/3];

switch PlotThis
    case 'FR'
        if strcmp(UnitFilt,'sparse')
            yminval = 0;
            ymaxval = 20;
        else
            yminval = 0;
            ymaxval = 40;
        end
    case 'FF'
        yminval = 0;
        ymaxval = 10;
    case 'TrV'
        yminval = 0;
        ymaxval = 0.7;
    case {'CV' 'xOWR'}
        yminval = 0;
        ymaxval = 4;
    case {'VS' 'Rel'}
        yminval = 0;
        ymaxval = 1;
    case 'MP'
        yminval = 0;
        ymaxval = 2*pi;
    case {'r_all' 'r_pdc' 'avTrCorr'}
        yminval = 0;
        ymaxval = 1;
    case 'covSpk'
        yminval = 0;
        ymaxval = 0.04;
    case {'CrSc' 'crW'}
        yminval = -0.2;
        ymaxval = 0.8;
    case 'OWR'
        if strcmp(UnitFilt,'sparse')
            yminval = -1.1;
            ymaxval = 8;
        else
            yminval = -1.1;
            ymaxval = 4;
        end
    otherwise
        keyboard
end

AMrates = [2 4 8 16 32];
ContextStr = {'Periodic' 'Irregular'};
SeqPosStr = {'Early' 'Late'};
DpTypeStr = {'Pdc-Pdc' 'IR-IR' 'Pdc-IR'};



%% Set up figures


% RATES SEPARATED
for ir = 1:numel(AMrates)
    
    hfr(ir) = figure;
    set(hfr(ir),'Position',fullscreen)
    for isp = 1:numel(DpTypeStr)
        hspr(ir,isp)=subplot(2,2*numel(DpTypeStr),[isp*2-1 isp*2]); hold on
        axis square
        plot([yminval ymaxval],[yminval ymaxval],':k','LineWidth',2);
        xlim([yminval ymaxval])
        ylim([yminval ymaxval])
        title(DpTypeStr{isp})
        ylabel([(PlotThis) ' of MPH 2'])
        xlabel([(PlotThis) ' of MPH 1'])
        hold off
    end
    hspr(ir,4)=subplot(2,2*numel(DpTypeStr),7:9); hold on
    plot([0 0],[0 40],':k','LineWidth',2);
    
    if any( strcmp( PlotThis, {'VS' 'MP' 'CrSc' 'crW' 'Rel'} ) ) 
        xlim([-ymaxval ymaxval])
    else
        xlim([-ymaxval ymaxval].*0.75)
    end
    
    ylim([0 40])
    xlabel([(PlotThis) ' MPH1 - ' (PlotThis) ' MPH2'])
    ylabel('Unit #')
    hold off
    
    hspr(ir,5)=subplot(2,2*numel(DpTypeStr),10:12); hold on
    plot([0 0],[0 40],':k','LineWidth',2);
    
    if any( strcmp( PlotThis, {'VS' 'MP' 'CrSc' 'crW' 'Rel'} ) )
        xlim([-ymaxval ymaxval])
    else
        xlim([-ymaxval ymaxval].*0.75)
    end
    
    ylim([0 1])
    xlabel([(PlotThis) ' MPH1 - ' (PlotThis) ' MPH2'])
    ylabel('Unit #')
    hold off
    
    suptitle(sprintf('%s of MPHs  |  stimuli: %s  |  %i Hz',PlotThis,select_stim,AMrates(ir)))
    
end
hold off



% OVERLAY

hfo = figure; hold on
set(hfo,'Position',fullscreen)

for isp = 1:numel(DpTypeStr)
    hspo(isp)=subplot(2,2*numel(DpTypeStr),[isp*2-1 isp*2]); hold on
    axis square
    plot([yminval ymaxval],[yminval ymaxval],':k','LineWidth',2);
    xlim([yminval ymaxval])
    ylim([yminval ymaxval])
    title(DpTypeStr{isp})
    ylabel([(PlotThis) ' of MPH 2'])
    xlabel([(PlotThis) ' of MPH 1'])
    hold off
end

hspo(4)=subplot(2,2*numel(DpTypeStr),7:9); hold on
plot([0 0],[0 40],':k','LineWidth',2);

if any(strcmp(PlotThis,{'VS' 'MP' 'CrSc' 'Rel'}))
    xlim([-ymaxval ymaxval])
else
    xlim([-ymaxval ymaxval].*0.75)
end

ylim([0 40])
xlabel([(PlotThis) ' MPH1 - ' (PlotThis) ' MPH2'])
ylabel('Unit #')

hspo(5)=subplot(2,2*numel(DpTypeStr),10:12); hold on
plot([0 0],[0 1],':k','LineWidth',2);

if any(strcmp(PlotThis,{'VS' 'MP' 'CrSc' 'Rel'}))
    xlim([-ymaxval ymaxval])
else
    xlim([-ymaxval ymaxval].*0.75)
end

ylim([0 0.4])
xlabel([(PlotThis) ' MPH1 - ' (PlotThis) ' MPH2'])
ylabel('Probability')

suptitle(sprintf('%s of MPHs  |  stimuli: %s  |  All AM rates',PlotThis,select_stim))
hold off




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











%%


% thisUnit = find(strcmp(UnitInfo.Session,'QB') & UnitInfo.Channel==3 & UnitInfo.Clu==32);


for iUn = 1:numel(UnitData)
    
    if isempty(UnitFilt) && strncmp(UnitInfo.RespType{iUn},'merged',2)
        continue
    elseif ~isempty(UnitFilt) && ~strncmp(UnitInfo.RespType{iUn},UnitFilt,2)
        continue
    end
    
    subject = UnitInfo(iUn,:).Subject{:};
    session = UnitInfo(iUn,:).Session{:}(1:2);
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
    
    % Get spiketimes and shift based on calculated integration time
    spiketimes = round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000);  %ms
    spiketimes = spiketimes-round(UnitData(iUn).IntTime_spk);
    
    % Get sound parameters
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials,dBSPL,LP,spiketimes);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    %%  IF THIS UNIT WAS SPLIT ACROSS SESSIONS, APPEND DATA TO MPH TABLE
    
    if numel(UnitData(iUn).Clu)>1
        
        session = UnitInfo(iUn,:).Session{:}(3:4);
        channel = UnitData(iUn).Channel(2);
        clu     = UnitData(iUn).Clu(2);
        
        % Load data files
        if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
            fprintf('Loading sess %s...\n',session)
            filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
            filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
        end
        if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
            filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
        end
        
        % Get spiketimes and shift based on calculated integration time
        spiketimes = round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000);  %ms
        spiketimes = spiketimes-round(UnitData(iUn).IntTime_spk);
        
        % Get sound parameters
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        MPH2 = makeMPHtable(TrialData,Info.artifact(channel).trials,dBSPL,LP,spiketimes);
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        % Merge the tables
        MPH = [MPH; MPH2];
        
    end
    
    
    
    
    
    
    
    
    
    
    %%  3) ANALYZE DATA
    
    ThisResp = UnitData(iUn);
    N = N+1;
    
    if strcmp(select_stim,'all')
        ThisResp.all = 1:5;
        %                     diffData_unit= [];
    end
    
    if ~isfield(ThisResp,select_stim)
        continue
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    anThisRate = AMrates(ThisResp.(select_stim));
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    for this_rate = anThisRate
        
        Period = 1000/this_rate;
        Theta  = linspace(0,2*pi,ceil(1000/this_rate));
        
        if isfield(ThisResp,'iBMF_FR') && ~isempty(ThisResp.iBMF_FR)
            logDistBMF = find(AMrates==this_rate)-ThisResp.iBMF_FR;
        else
            logDistBMF = nan;
        end
        
        
        %% ~~~~~~~~ Pdc-Pdc ~~~~~~~~                    where are the periods that follow 4 hz iti?
        
        subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate,:);
        PdStarttimes = unique(subMPH.Starttime)';
        
        % Skip Pdc periods that fall in first 100 ms
        if this_rate==8 || this_rate==16
            PdStarttimes(PdStarttimes<100) = [];
            
            % UNLESS it's 4 hz or 32 hz AND preceded by the same rate
        elseif this_rate==4  ||  this_rate==2  %or 2hz preceded by 4hz
            subMPH((subMPH.PrevStimID~=3) & (subMPH.Starttime<100),:)=[];
            PdStarttimes = unique(subMPH.Starttime)';
        elseif this_rate==32
            subMPH((subMPH.PrevStimID~=6) & (subMPH.Starttime<100),:)=[];
            PdStarttimes = unique(subMPH.Starttime)';
        end
        
        
        % Choose 2 or 4 random combinations of periods to
        % compare
        iPdPairs = randperm(length(PdStarttimes));
        PdPairs = PdStarttimes(iPdPairs(1:2));
        if this_rate>2
            PdPairs = [PdPairs; PdStarttimes(iPdPairs(3:4))];
            if this_rate>4
                PdPairs = [PdPairs; PdStarttimes(iPdPairs(5:6))];
            end
        end
        
        
        % Preallocate
        PdcStruct = struct;
        FR       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        SEM      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        FF       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        TrV      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        CV       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        VS       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        MP       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        CrSc     = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        Rel      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        covSpktimes = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        OWR      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        
        for ipair = 1:size(PdPairs,1)
            
            randPair = PdPairs(ipair,:);
            
            for ipd = 1:2
                
                for ipst = 1:numel(subMPH(subMPH.Starttime==randPair(ipd),:).raster)
                    
                    raster = subMPH(subMPH.Starttime==randPair(ipd),:).raster{ipst};
                    if ~(ceil(Period)==size(raster,2))
                        keyboard
                    end
                    
                    [CrSc(ipst),Rel(ipst)] = calculateCorrScore(raster,tVarBin,[num2str(this_rate) ' Hz period'],0);
                    
                    FR(1,ipst)  = sum(sum( raster ))/size(raster,1)/Period*1000;
                    SEM(1,ipst) = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                    FF(1,ipst)  = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                    TrV(1,ipst) = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                    CV(1,ipst)  = calc_binned_bootstrapped_CV(raster,tVarBin);
                    covSpktimes(ipst) = mean(diag(cov(raster)));
                    OWR(1,ipst) = ( mean(mean(raster(:,Theta<=OWbound),2)) - mean(mean(raster,2)) ) / mean(mean(raster,2));
                    
                    % Reformat raster for VS calculation
                    spktimes=[];
                    for it = 1:size(raster,1)
                        spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                    end
                    VS(1,ipst) = vectorstrength(spktimes,Period);
                    [MP(1,ipst),R] = meanphase(spktimes,Period);
                    if abs(R-VS(1,ipst))>1e-3, keyboard, end
                    
                end %ipst
                
                % Save data for this period
                
                PdcStruct(ipair,ipd).stimID = mode(subMPH.ThisStimID);
                PdcStruct(ipair,ipd).starttime = randPair(ipd);
                PdcStruct(ipair,ipd).FR  = mean(FR,'omitnan');
                PdcStruct(ipair,ipd).SEM = mean(SEM,'omitnan');
                PdcStruct(ipair,ipd).FF  = mean(FF,'omitnan');
                PdcStruct(ipair,ipd).TrV = mean(TrV,'omitnan');
                PdcStruct(ipair,ipd).CV  = mean(CV,'omitnan');
                PdcStruct(ipair,ipd).VS  = mean(VS,'omitnan');
                PdcStruct(ipair,ipd).MP  = mean(MP,'omitnan');
                PdcStruct(ipair,ipd).CrSc= mean(CrSc,'omitnan');
                PdcStruct(ipair,ipd).Rel = mean(Rel,'omitnan');
                PdcStruct(ipair,ipd).covSpk = mean(covSpktimes,'omitnan');
                PdcStruct(ipair,ipd).OWR = mean(OWR,'omitnan');
                
                PdcStruct(ipair,ipd).crW = 0;
                
            end %ipd
        end %ipair
        
        
        % Add Pdc data average for history comparison later
        
        %                     subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate & MPH.PrevPd==this_rate & MPH.Starttime>100,:);
        %                     if mode(subMPH.PrevPd)~=this_rate, keyboard, end
        
        meanPlotThis = mean([PdcStruct.(PlotThis)]);
        
        %                     HistoryData_Pdc = [HistoryData_Pdc; [ this_rate mode(subMPH.PrevPd) mean(subMPH.Prev500msFR) mean(subMPH.Prev100msFR) meanPlotThis ]];
        
        
        
        %% ~~~~~~~~  IR-IR  ~~~~~~~~
        
        
        subMPH = MPH(MPH.ThisStimID>6 & MPH.AMrate==this_rate,:);
        
        % Get all combinations of 2 periods (averaging across
        % prev stim)
        PdStarttimes = unique(subMPH.Starttime)';
        
        
        PdPairs = nchoosek(PdStarttimes,2);
        
        % Randomly remove a pair sometimes, to keep the sample
        % sizes similar across comparisons
        if numel(unique(subMPH.ThisStimID))>1 && rand(1,1)>0.333
            rnd_i_rm = randi(size(PdPairs,1));
            PdPairs(rnd_i_rm,:) = [];
        end
        
        
        % Preallocate
        IRStruct = struct;
        FR       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        SEM      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        FF       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        TrV      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        CV       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        VS       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        MP       = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        CrSc     = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        Rel      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        covSpktimes = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        OWR      = nan(1,numel(subMPH(subMPH.Starttime==PdStarttimes(1),:).raster));
        
        for ipair = 1:size(PdPairs,1)
            
            randPair = PdPairs(ipair,:);
            
            for ipd = 1:2
                
                for ipst = 1:numel(subMPH(subMPH.Starttime==randPair(ipd),:).raster)
                    
                    raster = subMPH(subMPH.Starttime==randPair(ipd),:).raster{ipst};
                    if ~(ceil(Period)==size(raster,2))
                        keyboard
                    end
                    
                    [CrSc(ipst),Rel(ipst)] = calculateCorrScore(raster,tVarBin,[num2str(this_rate) ' Hz period'],0);
                    
                    FR(1,ipst)  = sum(sum( raster ))/size(raster,1)/Period*1000;
                    SEM(1,ipst) = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                    FF(1,ipst)  = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                    TrV(1,ipst) = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                    CV(1,ipst)  = calc_binned_bootstrapped_CV(raster,tVarBin);
                    covSpktimes(ipst) = mean(diag(cov(raster)));
                    OWR(1,ipst) = ( mean(mean(raster(:,Theta<=OWbound),2)) - mean(mean(raster,2)) ) / mean(mean(raster,2));
                    
                    % Reformat raster for VS calculation
                    spktimes=[];
                    for it = 1:size(raster,1)
                        spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                    end
                    VS(1,ipst) = vectorstrength(spktimes,Period);
                    [MP(1,ipst),R] = meanphase(spktimes,Period);
                    if abs(R-VS(1,ipst))>1e-3, keyboard, end
                    
                end
                
                % Save data for this period
                
                IRStruct(ipair,ipd).stimID = mode(subMPH.ThisStimID);
                IRStruct(ipair,ipd).starttime = randPair(ipd);
                IRStruct(ipair,ipd).FR  = mean(FR,'omitnan');
                IRStruct(ipair,ipd).SEM = mean(SEM,'omitnan');
                IRStruct(ipair,ipd).FF  = mean(FF,'omitnan');
                IRStruct(ipair,ipd).TrV = mean(TrV,'omitnan');
                IRStruct(ipair,ipd).CV  = mean(CV,'omitnan');
                IRStruct(ipair,ipd).VS  = mean(VS,'omitnan');
                IRStruct(ipair,ipd).MP  = mean(MP,'omitnan');
                IRStruct(ipair,ipd).CrSc= mean(CrSc,'omitnan');
                IRStruct(ipair,ipd).Rel = mean(Rel,'omitnan');
                IRStruct(ipair,ipd).covSpk = mean(covSpktimes,'omitnan');
                IRStruct(ipair,ipd).OWR = mean(OWR,'omitnan');
                
                IRStruct(ipair,ipd).crW = 0;
                
                
            end %ipd
        end %ipair
        
        
        
        
        %% ~~~~~~ // Pdc-IR // ~~~~~~~~
        
        
        PdIRStruct = struct;
        theseIRs = unique(MPH.ThisStimID(MPH.ThisStimID>6))';
        
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
                MP       = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                CrSc     = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                Rel      = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                covSpktimes = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                OWR      = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                rAcrs    = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                rWith    = nan(1,numel(subMPH(subMPH.SeqPos==iseq,:).raster));
                
                for ipst = 1:numel(subMPH(subMPH.SeqPos==iseq,:).raster)
                    
                    
                    % Skip first period of AC if preceded by 32 Hz
                    if IRstim==7 && this_rate==32 && iseq==1 && ipst==5
                        continue
                    end
                    
                    % Skip first period of DB if preceded by 4 Hz
                    if IRstim==8 && this_rate==4 && iseq==1 && ipst==2
                        continue
                    end
                    
                    
                    raster = subMPH(subMPH.SeqPos==iseq,:).raster{ipst};
                    if ~(ceil(Period)==size(raster,2))
                        keyboard
                    end
                    
                    [CrSc(ipst),Rel(ipst)] = calculateCorrScore(raster,tVarBin,[num2str(this_rate) ' Hz period'],0);
                    
                    FR(ipst)  = sum(sum( raster ))/size(raster,1)/Period*1000;
                    SEM(ipst) = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                    FF(ipst)  = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                    TrV(ipst) = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                    CV(ipst)  = calc_binned_bootstrapped_CV(raster,tVarBin);
                    OWR(ipst) = ( mean(mean(raster(:,Theta<=OWbound),2)) - mean(mean(raster,2)) ) / mean(mean(raster,2));
%                     [r,p] = corrcoef(MPH_AllPds,sum(raster,1)/size(raster,1)*1000);
%                     rMPH_all(:,ipst) = [r(1,2); p(1,2)];
%                     [r,p] = corrcoef(MPH_PdcPds,sum(raster,1)/size(raster,1)*1000);
%                     rMPH_pdc(:,ipst) = [r(1,2); p(1,2)];
                    covSpktimes(ipst) = mean(diag(cov(raster)));
                    crW(ipst) = 0;
                    
                    % Reformat raster for VS calculation
                    spktimes=[];
                    for it = 1:size(raster,1)
                        spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                    end
                    VS(ipst) = vectorstrength(spktimes,Period);
                    [MP(ipst),R] = meanphase(spktimes,Period);
                    if abs(R-VS(ipst))>1e-3, keyboard, end
                    
                    % Add period to HistoryData
                    HistoryData = [HistoryData; ...
                        [ this_rate logDistBMF ...
                        subMPH(subMPH.SeqPos==iseq,:).PrevPd(ipst) ...
                        subMPH(subMPH.SeqPos==iseq,:).Prev500msFR(ipst) ...
                        subMPH(subMPH.SeqPos==iseq,:).Prev100msFR(ipst) ...
                        eval(sprintf('%s(ipst)',PlotThis))-meanPlotThis ] ];
                    
                    
                end %ipst
                
                % Average over previous stimuli and save for
                % comparison
                
                idx = 2*(IRstim-7) + iseq;
                
                PdIRStruct(idx,1).stimID = IRstim;
                PdIRStruct(idx,1).seqPos = iseq;
                PdIRStruct(idx,1).starttime = mode(subMPH(subMPH.SeqPos==iseq,:).Starttime);
                PdIRStruct(idx,1).FR   = mean(FR,'omitnan');
                PdIRStruct(idx,1).SEM  = mean(SEM,'omitnan');
                PdIRStruct(idx,1).FF   = mean(FF,'omitnan');
                PdIRStruct(idx,1).TrV  = mean(TrV,'omitnan');
                PdIRStruct(idx,1).CV   = mean(CV,'omitnan');
                PdIRStruct(idx,1).VS   = mean(VS,'omitnan');
                PdIRStruct(idx,1).MP   = mean(MP,'omitnan');
                PdIRStruct(idx,1).CrSc = mean(CrSc,'omitnan');
                PdIRStruct(idx,1).Rel = mean(Rel,'omitnan');
                PdIRStruct(idx,1).covSpk = mean(covSpktimes,'omitnan');
                PdIRStruct(idx,1).OWR = mean(OWR,'omitnan');
                
                
                % Correlation of MPHs within a context, from different previous
                % periods
                % Not sure this is the right comparison yet. First of all,
                % there are a vastly different number of cycles being
                % correlated across contexts and across rates. Probably would
                % be better to make a separate function for this analysis.
                [r,p] = corr(cell2mat(cellfun(@(x) mean(x,1), subMPH(subMPH.SeqPos==iseq,:).raster ,'UniformOutput',false))');
                PdIRStruct(idx,1).crW          = mean(r(~eye(size(r))));
                PdIRStruct(idx,1).corrWithin.r = mean(r(~eye(size(r))));
                PdIRStruct(idx,1).corrWithin.p = mean(p(~eye(size(p))));
                PdIRStruct(idx,1).corrWithin.n = numel(subMPH(subMPH.SeqPos==iseq,:).raster);
                
            end %iseq
            
            
        end %IRstim -- IRREGULAR first
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Now get matching PERIODIC context data
        
        subMPH = MPH(MPH.ThisStimID<7 & MPH.AMrate==this_rate,:);
        pdcstarttimes = unique(subMPH.Starttime)';
        
        % Skip Pdc periods that fall in first 100 ms
        if this_rate==2 || this_rate==8 || this_rate==16
            pdcstarttimes(pdcstarttimes<100) = [];
            
            % unless it's 4 hz or 32 hz preceded by the same rate
        elseif this_rate==4
            subMPH((subMPH.PrevStimID~=3) & (subMPH.Starttime<100),:)=[];
            pdcstarttimes = unique(subMPH.Starttime)';
        elseif this_rate==32
            subMPH((subMPH.PrevStimID~=6) & (subMPH.Starttime<100),:)=[];
            pdcstarttimes = unique(subMPH.Starttime)';
        end
        
        
        for iIRpd = 1:size(PdIRStruct,1)
            
            if isempty(PdIRStruct(iIRpd,1).stimID)
                continue
            end
            
            IRstim = PdIRStruct(iIRpd,1).stimID;
            iseq = PdIRStruct(iIRpd,1).seqPos;
            idx = 2*(IRstim-7) + iseq;
            if idx~=iIRpd, keyboard, end
            
            strtdiffs = PdIRStruct(idx,1).starttime - pdcstarttimes;
            %                         [~,ipd] = min(strtdiffs(strtdiffs>=0));
            [~,ipd] = min(abs(strtdiffs));
            
            % Preallocate
            FR      = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            SEM     = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            FF      = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            TrV     = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            CV      = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            VS      = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            MP      = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            CrSc    = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            Rel     = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            covSpktimes = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            OWR     = nan(1,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            %                         rMPH_all = nan(2,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            %                         rMPH_pdc = nan(2,numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster));
            
            for ipst = 1:numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster)
                
                raster = subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster{ipst};
                if ~(ceil(Period)==size(raster,2))
                    keyboard
                end
                
                [CrSc(ipst),Rel(ipst)] = calculateCorrScore(raster,tVarBin,[num2str(this_rate) ' Hz period'],0);
                
                FR(ipst)   = sum(sum( raster ))/size(raster,1)/Period*1000;
                SEM(ipst)  = std(sum(raster,2)/Period*1000) / sqrt(size(raster,1));
                FF(ipst)   = var(sum(raster,2)/Period*1000) / mean(sum(raster,2)/Period*1000);
                TrV(ipst)  = calc_binned_bootstrapped_TrTrVar(raster,tVarBin);
                CV(ipst)   = calc_binned_bootstrapped_CV(raster,tVarBin);
                covSpktimes(ipst) = mean(diag(cov(raster)));
                OWR(ipst)  = ( mean(mean(raster(:,Theta<=OWbound),2)) - mean(mean(raster,2)) ) / mean(mean(raster,2));
                %                             [r,p] = corrcoef(MPH_AllPds,sum(raster,1)/size(raster,1)*1000);
                %                             rMPH_all(:,ipst) = [r(1,2); p(1,2)];
                %                             [r,p] = corrcoef(MPH_PdcPds,sum(raster,1)/size(raster,1)*1000);
                %                             rMPH_pdc(:,ipst) = [r(1,2); p(1,2)];
                
                % Reformat raster for VS calculation
                spktimes=[];
                for it = 1:size(raster,1)
                    spktimes = [spktimes find(raster(it,:)) + (it-1)*Period];
                end
                VS(ipst) = vectorstrength(spktimes,Period);
                [MP(ipst),R] = meanphase(spktimes,Period);
                if abs(R-VS(ipst))>1e-3, keyboard, end
                
            end %ipst
            
            
            PdIRStruct(idx,2).stimID = this_rate;
            PdIRStruct(idx,2).seqPos = iseq;
            PdIRStruct(idx,2).starttime = mode(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).Starttime);
            PdIRStruct(idx,2).FR  = mean(FR,'omitnan');
            PdIRStruct(idx,2).SEM = mean(SEM,'omitnan');
            PdIRStruct(idx,2).FF  = mean(FF,'omitnan');
            PdIRStruct(idx,2).TrV = mean(TrV,'omitnan');
            PdIRStruct(idx,2).CV  = mean(CV,'omitnan');
            PdIRStruct(idx,2).VS  = mean(VS,'omitnan');
            PdIRStruct(idx,2).MP  = mean(MP,'omitnan');
            PdIRStruct(idx,2).CrSc= mean(CrSc,'omitnan');
            PdIRStruct(idx,2).Rel = mean(Rel,'omitnan');
            PdIRStruct(idx,2).covSpk = mean(covSpktimes,'omitnan');
            PdIRStruct(idx,2).OWR = mean(OWR,'omitnan');
            
            % Correlation of MPHs within a context, from different previous
            % periods
            % Not sure this is the right comparison yet. First of all,
            % there are a vastly different number of cycles being
            % correlated across contexts and across rates. Probably would
            % be better to make a separate function for this analysis.
            % Also, should be binned at lower resolution than 1 ms.
            [r,p] = corr(cell2mat(cellfun(@(x) mean(x,1), subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster ,'UniformOutput',false))');
            PdIRStruct(idx,2).crW          = mean(r(~eye(size(r))));
            PdIRStruct(idx,2).corrWithin.r = mean(r(~eye(size(r))));
            PdIRStruct(idx,2).corrWithin.p = mean(p(~eye(size(p))));
            PdIRStruct(idx,2).corrWithin.n = numel(subMPH(subMPH.Starttime==pdcstarttimes(ipd),:).raster);
            
            
            
        end %iIRpd
        
        
        
        %% ========  PLOT  =======
        
        
        % RATES SEPARATED
        
        figure(hfr(AMrates==this_rate)); hold on
        
        % Pdc-Pdc
        subplot(hspr(AMrates==this_rate,1)); hold on
        %                     if this_rate<8
        scatter([PdcStruct(:,1).(PlotThis)],[PdcStruct(:,2).(PlotThis)],dotsize,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
        %                     else
        %                         MPH1range = prctile([PdcData(:,1).(PlotThis)],[5 95]);
        %                         MPH2range = prctile([PdcData(:,2).(PlotThis)],[5 95]);
        %                         plot(MPH1range,[median([PdcData(:,2).(PlotThis)]) median([PdcData(:,2).(PlotThis)])],'-b')
        %                         plot([median([PdcData(:,1).(PlotThis)]) median([PdcData(:,1).(PlotThis)])],MPH2range,'-b')
        %                     end
        hold off
        
        % IR-IR
        subplot(hspr(AMrates==this_rate,2)); hold on
        scatter([IRStruct(:,1).(PlotThis)],[IRStruct(:,2).(PlotThis)],dotsize,'o','MarkerFaceColor','g','MarkerEdgeColor','g')
        hold off
        
        % Pdc-IR
        subplot(hspr(AMrates==this_rate,3)); hold on
        scatter([PdIRStruct(:,2).(PlotThis)], [PdIRStruct(:,1).(PlotThis)], dotsize,'o',...
            'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:))%,...
        %                         'MarkerFaceAlpha',alphval,'MarkerEdgeAlpha',alphval)
        hold off
        
        % Overlay -- difference plot
        subplot(hspr(AMrates==this_rate,4)); hold on
        %                     if this_rate<8
        scatter([PdcStruct(:,1).(PlotThis)]-[PdcStruct(:,2).(PlotThis)], N*ones(size([PdcStruct(:,1).(PlotThis)])),...
            dotsize,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
        %                     else
        %                         diffPctl = prctile([PdcData(:,1).(PlotThis)]-[PdcData(:,2).(PlotThis)],[5 95]); %for hist of diffs
        %                         plot(diffPctl, [N N], 'b','LineWidth',3)
        %                     end
        
        scatter([PdIRStruct(:,2).(PlotThis)]-[PdIRStruct(:,1).(PlotThis)], N*ones(size([PdIRStruct(:,1).(PlotThis)])), dotsize,'o',...
            'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',[1 1 1])
        ylim([0 N+1])
        hold off
        
        
        
        % OVERLAY
        
        figure(hfo); hold on
        
        % Pdc-Pdc
        subplot(hspo(1)); hold on
        scatter([PdcStruct(:,1).(PlotThis)],[PdcStruct(:,2).(PlotThis)],dotsize,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
        hold off
        
        % IR-IR
        subplot(hspo(2)); hold on
        scatter([IRStruct(:,1).(PlotThis)],[IRStruct(:,2).(PlotThis)],dotsize,'o','MarkerFaceColor','g','MarkerEdgeColor','g')
        hold off
        
        % Pdc-IR
        subplot(hspo(3)); hold on
        scatter([PdIRStruct(:,2).(PlotThis)], [PdIRStruct(:,1).(PlotThis)], dotsize,'o',...
            'MarkerFaceColor',colors(subMPH.ThisStimID(1),:),'MarkerEdgeColor',colors(subMPH.ThisStimID(1),:))%,...
        hold off
        
        
        
        %% ------ Store data  ------
        
        % All data collected
        GroupData(N).(['rate' num2str(this_rate)]).PdcData   = PdcStruct;
        GroupData(N).(['rate' num2str(this_rate)]).IRData    = IRStruct;
        GroupData(N).(['rate' num2str(this_rate)]).PdIRData  = PdIRStruct;
        
        
        % For stats of Pdc-IR, all AM rates overlayed
        Pdcmat    = [ Pdcmat;  [ [PdcStruct(:,1).(PlotThis)]'  [PdcStruct(:,2).(PlotThis)]' ] ];
        IRmat     = [ IRmat;   [ [IRStruct(:,1).(PlotThis)]'   [IRStruct(:,2).(PlotThis)]'  ] ];
        PdIRmat   = [ PdIRmat; [ [PdIRStruct(:,2).(PlotThis)]' [PdIRStruct(:,1).(PlotThis)]'] ];
        
        
        % Pdc-IR differences, as a function of stimulus
        if strcmp(select_stim,'all')
            
            % By AM rate
            hfdrData = [hfdrData; find(AMrates==this_rate) PdIRStruct(idx,1).(PlotThis)-PdIRStruct(idx,2).(PlotThis)];
            
            % By distance to BMF-fr
            if isfield(ThisResp,'iBMF_FR') && ~isempty(ThisResp.iBMF_FR)
                hfddData = [hfddData; logDistBMF PdIRStruct(idx,1).(PlotThis)-PdIRStruct(idx,2).(PlotThis)];
            end
        end
        
        
        
        
    end %this_rate
    
    
    
    %% Plot overlay now, with Pdc-Pdc datapoints on bottom layer
    
    % Overlay -- difference plot
    figure(hfo); hold on
    subplot(hspo(4)); hold on
    
    for ir=anThisRate
        scatter([GroupData(N).(['rate' num2str(ir)]).PdcData(:,1).(PlotThis)]-[GroupData(N).(['rate' num2str(ir)]).PdcData(:,2).(PlotThis)],...
            N*ones(size([GroupData(N).(['rate' num2str(ir)]).PdcData(:,1).(PlotThis)])),...
            dotsize/1.5,'o','MarkerFaceColor','b','MarkerEdgeColor','b')
    end
    for ir=anThisRate
        scatter([GroupData(N).(['rate' num2str(ir)]).PdIRData(:,2).(PlotThis)]-[GroupData(N).(['rate' num2str(ir)]).PdIRData(:,1).(PlotThis)],...
            N*ones(size([GroupData(N).(['rate' num2str(ir)]).PdIRData(:,1).(PlotThis)])),...
            dotsize/3,'o','MarkerFaceColor',colors(1+find(ir==AMrates),:),'MarkerEdgeColor',[1 1 1])
    end
    
    ylim([0 N+1])
    hold off
    
    
    
    
    
    %%   End of datapoint
    
end %iUn


%% Histogram of differences for Pdc-Pdc and Pdc-IR comparisons, and stats


% For overlay of all rates first
PdcDiffs  = Pdcmat(:,1) - Pdcmat(:,2);
PdIRDiffs = PdIRmat(:,1) - PdIRmat(:,2);

p_rs=ranksum(PdcDiffs,PdIRDiffs);

h_l(1) = lillietest(PdcDiffs);
h_l(2) = lillietest(PdIRDiffs);
if any(h_l==0)
    [~,p_f] = deal(nan);
else
    [~,p_f] = vartest2(PdcDiffs,PdIRDiffs);
end


histbins = -ymaxval:(ymaxval/20):ymaxval;

figure(hfo); hold on
subplot(hspo(5)); hold on

histogram(PdcDiffs,histbins,'FaceColor','b','EdgeColor','b','Normalization','probability')
hold on
% histogram(IRDiffs,histbins,'FaceColor','g','EdgeColor','g','Normalization','probability')
histogram(PdIRDiffs,histbins,'FaceColor','r','EdgeColor','r','Normalization','probability')
title(sprintf('ranksum p = %0.2e,  F-test p=%0.2e',p_rs,p_f))
hold off


% Now for each AM rate separately
for ir = 1:5
    
    PdcDiffs=[];
    PdIRDiffs=[];

    for in = 1:numel(GroupData)
        if ~isfield(GroupData(in),['rate' num2str(AMrates(ir))]) || isempty(GroupData(in).(['rate' num2str(AMrates(ir))]))
            continue
        end
        % Pdc - Pdc
        PdcDiffs  = [PdcDiffs;  [GroupData(in).(['rate' num2str(AMrates(ir))]).PdcData(:,2).(PlotThis)]' - [GroupData(in).(['rate' num2str(AMrates(ir))]).PdcData(:,1).(PlotThis)]' ];
        % Pdc - IR
        PdIRDiffs = [PdIRDiffs; [GroupData(in).(['rate' num2str(AMrates(ir))]).PdIRData(:,2).(PlotThis)]' - [GroupData(in).(['rate' num2str(AMrates(ir))]).PdIRData(:,1).(PlotThis)]' ];
    end 
    
    if isempty(PdcDiffs) || isempty(PdIRDiffs)
        continue
    end
    
    % Stats
    p_rs=ranksum(PdcDiffs,PdIRDiffs);
    
    try
        h_l(1) = lillietest(PdcDiffs);
    catch
        h_l(1) = nan;
    end
    try
        h_l(2) = lillietest(PdIRDiffs);
    catch 
        h_l(2) = nan;
    end
    if any(h_l==0)
        [~,p_f] = deal(nan);
    else
        [~,p_f] = vartest2(PdcDiffs,PdIRDiffs);
    end
    
    % Plot
    figure(hfr(ir));
    subplot(hspr(ir,5)); hold on
    
    histogram(PdcDiffs,histbins,'FaceColor','b','EdgeColor','b','Normalization','probability')
    hold on
    histogram(PdIRDiffs,histbins,'FaceColor','r','EdgeColor','r','Normalization','probability')
    title(sprintf('ranksum p = %0.2e,  F-test p=%0.2e',p_rs,p_f))
    ylim([0 0.4])
    hold off
end



%% Stats for Pdc-IR overlay plots


% Pdc-Pdc
idx = ~isnan(Pdcmat(:,2).*Pdcmat(:,1));
[r,p_r] = corrcoef( Pdcmat(idx,1), Pdcmat(idx,2) );
p_wsr   = signrank( Pdcmat(:,1), Pdcmat(:,2) );

figure(hfo); hold on
subplot(hspo(1)); hold on
title(sprintf('Pdc-Pdc comparison\nWSR p=%0.2e, Pearson r=%0.2f, p=%0.2e',p_wsr,r(1,2),p_r(1,2)))
hold off

% IR-IR 
idx = ~isnan(IRmat(:,2).*IRmat(:,1));
[r,p_r] = corrcoef( IRmat(idx,2), IRmat(idx,1) );
p_wsr = signrank( IRmat(:,2), IRmat(:,1) );

subplot(hspo(2)); hold on
title(sprintf('IR-IR comparison\nWSR p=%0.2e, Pearson r=%0.2f, p=%0.2e',p_wsr,r(1,2),p_r(1,2)))
hold off

% Pdc-IR 
idx = ~isnan(PdIRmat(:,2).*PdIRmat(:,1));
[r,p_r] = corrcoef( PdIRmat(idx,2), PdIRmat(idx,1) );
p_wsr = signrank( PdIRmat(:,2), PdIRmat(:,1) );

subplot(hspo(3)); hold on
title(sprintf('Pdc-IR comparison\nWSR p=%0.2e, Pearson r=%0.2f, p=%0.2e',p_wsr,r(1,2),p_r(1,2)))
hold off




%% IR-Pdc diff against AM rate/distance to BMF
if strcmp(select_stim,'all')
    
    % Diff IR-Pdc, as a function of AMrate
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
    hfdd = figure; hold on
    set(hfdd,'Position',narrow)
    plot([-5 5],[0 0],':k');
    xlim([-5 5])
    ylim((ymaxval)*[-1 1])
    ylabel('diff IR-Pdc')
    xlabel('log distance to BMF')
    hold off
    
    
    figure(hfdr); hold on
    plotSpread(hfdrData(:,2),'distributionIdx',hfdrData(:,1),'distributionColors','k')
    boxplot(hfdrData(:,2)',hfdrData(:,1)','Positions',unique(hfdrData(:,1))',...
        'Colors','r','notch','on','Symbol','','Labels',Info.stim_ID_key(2:6)')
    p = kruskalwallis(hfdrData(:,2)',hfdrData(:,1)','off');
    title(sprintf('%s context diff as a fct of AM rate\nKW p=%0.2e',PlotThis,p))
    hold off
    
    figure(hfdd); hold on
    
    plotSpread(hfddData(:,2),'distributionIdx',hfddData(:,1),'distributionColors','k')
    boxplot(hfddData(:,2)',hfddData(:,1)','Positions',unique(hfddData(:,1))','Colors','r','notch','on','Symbol','')
    p = kruskalwallis(hfddData(:,2)',hfddData(:,1)','off');
    title(sprintf('%s context diff as a fct of distance from BMF-fr\nKW p=%0.2e',PlotThis,p))
    hold off
end



%% Check for a relationship between IR/Pdc response and prev pd and prev FR

% HistoryData   = [ AMrate  PrevPd  PrevFR500  PrevFR100  DependentVar_diff_from_Pdc ]

if strcmp(UnitFilt,'sparse')
    setxlim = [0 15];
else
    setxlim = [0 40];
end

hfd(1)=figure;
set(gcf,'Position',widerect)
hold on
for ir=1:5 
    subplot(1,5,ir); 
    plot(setxlim,[0 0],'k'); hold on
    
    ihf  = HistoryData(:,1)==AMrates(ir);
    
    plot( HistoryData(ihf,4), HistoryData(ihf,6), 'ok','LineWidth',2,'MarkerSize',10)
    
    title([num2str(AMrates(ir)) ' Hz period'])
    set(gca,'xlim',setxlim,'xscale','linear','ylim',ymaxval*[-1 1])
    hold off
end
suptitle('prev 500 ms FR')
hold off


hfd(2)=figure;
set(gcf,'Position',widerect)
hold on
for ir=1:5 
    subplot(1,5,ir); hold on
    plot(setxlim,[0 0],'k')
    
    ihf  = HistoryData(:,1)==AMrates(ir);
    
    plot( HistoryData(ihf,5), HistoryData(ihf,6), 'ok','LineWidth',2,'MarkerSize',10)
    
    title([num2str(AMrates(ir)) ' Hz period'])
    set(gca,'xlim',setxlim,'xscale','linear','ylim',ymaxval*[-1 1])
    hold off
end
suptitle('prev 100 ms FR')
hold off


setxlim = [1 64];

hfd(3)=figure;
set(gcf,'Position',widerect)
hold on
for ir=1:5 
    subplot(1,5,ir); hold on
    plot(setxlim,[0 0],'k')

    ihf  = HistoryData(:,1)==AMrates(ir);
    
    plot( HistoryData(ihf,3), HistoryData(ihf,6), 'ok','LineWidth',2,'MarkerSize',10)
    
    title([num2str(AMrates(ir)) ' Hz period'])
    set(gca,'xlim',setxlim,'xscale','log','ylim',ymaxval*[-1 1])
    hold off
end
suptitle('prev Pd, effect resp to AM rate')
hold off



dBMF = unique(HistoryData(~isnan(HistoryData(:,2)),2));

hfd(4)=figure;
set(gcf,'Position',widerect)
hold on
for ir=1:length(dBMF)
    subplot(1,length(dBMF),ir); hold on
    plot(setxlim,[0 0],'k')

    ihf  = HistoryData(:,2)==dBMF(ir);
    
    plot( HistoryData(ihf,3), HistoryData(ihf,6), 'ok','LineWidth',2,'MarkerSize',10)
    
    title(['BMF + ' num2str(dBMF(ir)) ])
    set(gca,'xlim',setxlim,'xscale','log','ylim',ymaxval*[-1 1])
    hold off
end
suptitle('prev Pd, effect on resp by dBMF')
hold off


setxlim = [min(dBMF)-1 max(dBMF)+1];

hfd(5)=figure;
set(gcf,'Position',widerect)
hold on
for ir=1:5 
    subplot(1,5,ir); hold on
    plot(setxlim,[0 0],'k')

    ihf  = (HistoryData(:,1)==AMrates(ir)) & ~isnan(HistoryData(:,2));
    
    plot( HistoryData(ihf,2), HistoryData(ihf,6), 'ok','LineWidth',2,'MarkerSize',10)
    
    title([num2str(AMrates(ir)) ' Hz period'])
    set(gca,'xlim',setxlim,'ylim',ymaxval*[-1 1])
    hold off
end
suptitle('This Pd relation to BMF')
hold off








%% SAVE PLOTS AND DATA

if isempty(UnitFilt)
    UnitFilt = 'All';
end

savedir = fullfile(fn.processed,'MPHcompare',UnitFilt);
if ~exist(savedir,'dir')
    mkdir([savedir '/svg']);
    mkdir([savedir '/eps']);
end

for ir = 1:5
    savename = sprintf('MPHs_%iHz_%s_%s',AMrates(ir),PlotThis,select_stim);
    % eps
    print_eps_kp(hfr(ir),fullfile([savedir '/eps'],[savename '_rates']))
    % svg
    print_svg_kp(hfr(ir),fullfile([savedir '/svg'],[savename '_rates']))
end

savename = sprintf('MPHs_ALLrates_%s_%s',PlotThis,select_stim);
print_eps_kp(hfo,fullfile([savedir '/eps'],[savename '_overlay']))
print_svg_kp(hfo,fullfile([savedir '/svg'],[savename '_overlay']))


% Also save comparison plots

savedir = [savedir '_Comparisons'];
if ~exist(savedir,'dir')
    mkdir([savedir '/svg']);
    mkdir([savedir '/eps']);
end

if strcmp(select_stim,'all')
    savename = sprintf('DiffByAMrate_%s_%s',PlotThis,select_stim);
    print_eps_kp(hfdr,fullfile([savedir '/eps'],savename))
    print_svg_kp(hfdr,fullfile([savedir '/svg'],savename))
    
    savename = sprintf('DiffByBMF_%s_%s',PlotThis,select_stim);
    print_eps_kp(hfdd,fullfile([savedir '/eps'],savename))
    print_svg_kp(hfdd,fullfile([savedir '/svg'],savename))
end

savename = sprintf('PdDiffs_FR500_%s_%s',PlotThis,select_stim);
print_eps_kp(hfd(1),fullfile([savedir '/eps'],savename))
print_svg_kp(hfd(1),fullfile([savedir '/svg'],savename))

savename = sprintf('PdDiffs_FR100_%s_%s',PlotThis,select_stim);
print_eps_kp(hfd(2),fullfile([savedir '/eps'],savename))
print_svg_kp(hfd(2),fullfile([savedir '/svg'],savename))

savename = sprintf('PdDiffs_PrevPd_%s_%s',PlotThis,select_stim);
print_eps_kp(hfd(3),fullfile([savedir '/eps'],savename))
print_svg_kp(hfd(3),fullfile([savedir '/svg'],savename))

savename = sprintf('PdDiffs_dBMF_PrevPd_%s_%s',PlotThis,select_stim);
print_eps_kp(hfd(4),fullfile([savedir '/eps'],savename))
print_svg_kp(hfd(4),fullfile([savedir '/svg'],savename))

savename = sprintf('PdDiffs_ThisPd_dBMF_%s_%s',PlotThis,select_stim);
print_eps_kp(hfd(5),fullfile([savedir '/eps'],savename))
print_svg_kp(hfd(5),fullfile([savedir '/svg'],savename))




% Also save GroupData
savename = sprintf('MPHs_%s',select_stim);
save(fullfile(savedir,savename),'GroupData','-v7.3')





% keyboard





end %function






