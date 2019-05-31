function tuningCurveComparisons
%
%  tuningCurveComparisons
%
%
%
%  KP, 2019-04
%


global fn AMrates rateVec_AC rateVec_DB trMin RateStream
fn = set_paths_directories([],[],1);

%!!!!!!!!!!!!!!!!!
trMin   =  10;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!
N=0;


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

AMrates = [2 4 8 16 32];


USE_MEASURE = 'FR';
switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 0.01;
                axmax = 10*ceil(max([UnitData.BaseFR])/10);
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
colors = [ 0 200 150;...
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


%% Select a datapoint (or a few?) to highlight

ex_subj = 'xWWWf_253400';
ex_sess = 'GA';
ex_ch   = 16;
ex_clu  = 1;


%% Preallocate

N=0;
NIR=0;

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
    
    close all
    
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
%     if UnitData(iUn).kw_p>0.05
%         continue
%     end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    subjcol     = [1 1 1];
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
    
    spkshift    = UnitData(iUn).IntTime_spk;
    
    
    % Load data files
    
%     if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
%     end
%     if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
%     end
    if ~isfield(Info,'artifact')
        continue
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
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    %%  3) PLOT DATA
    
    T = 250;
    
    Pdc_Nspikes_mean = nan(1,5);
    Pdc_Nspikes_sem  = nan(1,5);
    Pdc_FR_mean      = nan(1,5);
    Pdc_FR_sem       = nan(1,5);
    
    for this_rate = AMrates
        
        Period = 1000/this_rate;
        irate     = find(AMrates==this_rate);
        MPH_rate  = MPH(MPH.AMrate==this_rate,:);
        
        % Periodic indices 
        iPdc = find( (MPH_rate.ThisStimID==irate+1 & MPH_rate.Starttime>=T)...
            | (MPH_rate.ThisStimID==irate+1 & MPH_rate.PrevStimID==irate+1) );
        if ~all(MPH_rate.PrevPd(iPdc)==this_rate), keyboard, end
        
        this_MPH_raster     = [];
        for ip = iPdc'
            this_MPH_raster = [this_MPH_raster; MPH_rate(ip,:).raster{:}];
        end 
        
        Pdc_Nspikes_mean(irate) = mean(sum(this_MPH_raster,2));
        Pdc_Nspikes_sem(irate)  = std(sum(this_MPH_raster,2)) / sqrt(size(this_MPH_raster,1));
        
        Pdc_FR_mean(irate)      = mean( sum( this_MPH_raster,2 )/Period*1000 );
        Pdc_FR_sem(irate)       = std(  sum( this_MPH_raster,2 )/Period*1000 ) / sqrt(size(this_MPH_raster,1));
    end
    
    % Plot periodic MTF
    hf   =   figure;
    set(hf,'Position',fullscreen)
    errorbar(Pdc_FR_mean,Pdc_FR_sem,'k','LineWidth',3)
    xlim([0 6])
    hold on
    
    
        
    for this_rate = AMrates
        
        % Set up
        Period = 1000/this_rate;
        PdData = struct;
        T = 250;
        
        irate = find(AMrates==this_rate);
        MPH_rate = MPH(MPH.AMrate==this_rate,:);
        
        
        % Irregular context info
        theseIRs = unique(MPH_rate.ThisStimID(MPH_rate.ThisStimID>6));
        iIR = []; nIR = numel(theseIRs); 
        np=0;
        for thisIR = theseIRs
            iIR = [iIR find(MPH_rate.Starttime>=T & MPH_rate.ThisStimID==thisIR) ];
            PPds = unique(MPH_rate.PrevPd( MPH_rate.Starttime>=T & MPH_rate.ThisStimID==thisIR ));
            for ip = PPds'
                np = np+1;
            end
        end
        
        
        
        %==================================================================
        %% Now the pds that come right after a transition
        
        iTrans = find(MPH_rate.ThisStimID==irate+1 & MPH_rate.Starttime<T/2 & MPH_rate.PrevStimID~=irate+1 )';
        PStIDs = [MPH_rate(iTrans,:).PrevStimID]';
        
        for ips = unique(PStIDs)
            
            this_MPH_raster = [];
            this_MPH_FR     = [];
            np=np+1;
            
            for ip = iTrans(PStIDs==ips)
                
                this_MPH_raster = [this_MPH_raster; MPH_rate(ip,:).raster{:}];
                this_MPH_FR     = [this_MPH_FR;     MPH_rate(ip,:).FRsmooth{:}];
                
            end %iPdc
            
            if ~(ceil(Period)==size(this_MPH_raster,2))
                keyboard
            end
            
            % Reformat raster for VS calculation
            spktimes=[];
            for it = 1:size(this_MPH_raster,1)
                spktimes = [spktimes find(this_MPH_raster(it,:))];
            end
            VS = vectorstrength(spktimes,Period);
            
            % Store data
            PdData(np,1).stimID = MPH_rate(ip,:).ThisStimID;
            PdData(np,1).FR  = sum(sum( this_MPH_raster ))/size(this_MPH_raster,1)/Period*1000;
            PdData(np,1).SEM = std(sum(this_MPH_raster,2)/Period*1000) / sqrt(size(this_MPH_raster,1));
            PdData(np,1).FF  = var(sum(this_MPH_raster,2)/Period*1000) / mean(sum(this_MPH_raster,2)/Period*1000);
            PdData(np,1).TrV = calc_binned_bootstrapped_TrTrVar(this_MPH_raster,tVarBin);
            PdData(np,1).VS  = VS;
            
            % Add to plot
            plot(irate, mean(sum(this_MPH_raster,2)/Period*1000 ),'o','MarkerSize',15,'LineWidth',2,...
                'Color',colors(ips,:))
            
        end %PStIDs
        
        
        %==================================================================
        %% Get IRREGULAR context data
        
        for thisIR = theseIRs
            
            PPds = unique(MPH_rate.PrevPd( MPH_rate.Starttime>T & MPH_rate.ThisStimID==thisIR ));
            
            for ipp = PPds'
                
                np = np+1;
                
                this_MPH_raster = [];
                this_MPH_FR     = [];
                iIR = find(MPH_rate.PrevPd==ipp & MPH_rate.Starttime>T & MPH_rate.ThisStimID==thisIR)';
                
                for ii = iIR
                    this_MPH_raster = [this_MPH_raster; MPH_rate(ii,:).raster{:}];
                    this_MPH_FR     = [this_MPH_FR;     MPH_rate(ii,:).FRsmooth{:}];
                end
                
                % Reformat raster for VS calculation
                spktimes=[];
                for it = 1:size(this_MPH_raster,1)
                    spktimes = [spktimes find(this_MPH_raster(it,:))];
                end
                VS = vectorstrength(spktimes,Period);
                
                % Store data
                PdData(np,1).stimID = MPH_rate(ii,:).ThisStimID;
                PdData(np,1).FR  = sum(sum( this_MPH_raster ))/size(this_MPH_raster,1)/Period*1000;
                PdData(np,1).SEM = std(sum(this_MPH_raster,2)/Period*1000) / sqrt(size(this_MPH_raster,1));
                PdData(np,1).FF  = var(sum(this_MPH_raster,2)/Period*1000) / mean(sum(this_MPH_raster,2)/Period*1000);
                PdData(np,1).TrV = calc_binned_bootstrapped_TrTrVar(this_MPH_raster,tVarBin);
                PdData(np,1).VS  = VS;
                
                
                % Add to plot
                plot(irate,mean(sum(this_MPH_raster,2)/Period*1000 ),'x','MarkerSize',15,'LineWidth',2,...
                    'Color',colors(1+find(ipp==AMrates),:))
                
            end %PPds
        end %thisIR
        
    end  %this_rate
    
end %iUn

keyboard





end %function




