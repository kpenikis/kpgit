function plotMPH_aggregate(SUBJECT,SESSION)
%
%  plotMPH_aggregate(SUBJECT,SESSION)
%   Plots MPHs from aggregated spiketimes of RS and NS units. 
%   *** So far, only tested with Units file from single session !!
%
%  KP, 2019-03
%


close  all

global fn AMrates rateVec_AC rateVec_DB trMin RateStream
fn = set_paths_directories([],[],1);

%!!!!!!!!!!!!!!!!!
trMin   =  10;
%!!!!!!!!!!!!!!!!!
tVarBin = 31;
%!!!!!!!!!!!!!!!!!
yminval = 0;
ymaxval = 100;
%!!!!!!!!!!!!!!!!!
N=0;


%% Load data

[UnitInfo, UnitData, Info, TrialData, Clusters, ~, artifactTrs ] = collectRasterDataSession(SUBJECT,SESSION);
[spiketimes_NS, spiketimes_RS] = aggregateNSRSspikes(UnitData,UnitInfo,Clusters);


%% Set up figure

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
halfscreen = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


AMrates = [2 4 8 16 32];
ContextStr = {'Periodic' 'IR (AC)' 'IR (DB)'};
SeqPosStr = {'Early' 'Late'};

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

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


%% Now get MPH data for each type

for UnType = {'NS' 'RS'}
    
    switch UnType{:}
        case 'NS'
            spiketimes = spiketimes_NS;
        case 'RS'
            spiketimes = spiketimes_RS;
    end
    
    
    %%
    % Get sound parameters
    [dBSPL,LP] = theseSoundParams(TrialData);
    if numel(dBSPL)>1 || numel(LP)>1
        keyboard
    end
    
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,artifactTrs,dBSPL,LP,spiketimes);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    %%  3) PLOT DATA
    
    for this_rate = AMrates
        
        % Set up
        Period = 1000/this_rate;
        PdData = struct;
        T = 250;
        
        irate = find(AMrates==this_rate);
        MPH_rate = MPH(MPH.AMrate==this_rate,:);
           
        % Periodic indices first
        iPdc = find( (MPH_rate.ThisStimID==irate+1 & MPH_rate.Starttime>=T)...
            | (MPH_rate.ThisStimID==irate+1 & MPH_rate.PrevStimID==irate+1) );
        if ~all(MPH_rate.PrevPd(iPdc)==this_rate), keyboard, end
        np=1;
        
        % Irregular context info
        theseIRs = unique(MPH_rate.ThisStimID(MPH_rate.ThisStimID>6));
        iIR = []; nIR = numel(theseIRs);
        for thisIR = theseIRs
            iIR = [iIR find(MPH_rate.Starttime>=T & MPH_rate.ThisStimID==thisIR) ];
            PPds = unique(MPH_rate.PrevPd( MPH_rate.Starttime>=T & MPH_rate.ThisStimID==thisIR ));
            for ip = PPds'
                np = np+1;
            end
        end
        
        
        %==========================================================
        %% Set up figure
        %==========================================================
        
        
        hfe(irate)   =   figure;
        set(hfe(irate),'Position',fullscreen)
        
        for ir=1:2
            for ic = 1:3
                
                hsp(ir,ic)=subplot(2,3,ir*3-3+ic); hold on
                axis square
                
                if ir==2
                    xlabel('Time in period')
                end
                if ic==1
                    ylabel('Spikes/sec')
                end
            end
        end
        linkaxes(hsp,'xy')
        xlim([0 ceil(1000/this_rate)])
        ylim([yminval ymaxval])
        
        
        %==================================================================
        %% Get PERIODIC context first
        
        this_MPH_raster     = [];
        this_MPH_FR         = [];
        np = 1;
        
        for ip = iPdc'
            
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
        PdData(1,1).stimID = MPH_rate(ip,:).ThisStimID;
        PdData(1,1).FR  = sum(sum( this_MPH_raster ))/size(this_MPH_raster,1)/Period*1000;
        PdData(1,1).SEM = std(sum(this_MPH_raster,2)/Period*1000) / sqrt(size(this_MPH_raster,1));
        PdData(1,1).FF  = var(sum(this_MPH_raster,2)/Period*1000) / mean(sum(this_MPH_raster,2)/Period*1000);
        PdData(1,1).TrV = calc_binned_bootstrapped_TrTrVar(this_MPH_raster,tVarBin);
        PdData(1,1).VS  = VS;
        
        
        % - - - - - - - - - - - - - - - - 
        %           PLOT Pdc
        
        histbin   = floor(size(this_MPH_raster,2)/30);
        histMP    = []; histMP = binspikecounts(mean(this_MPH_raster,1),histbin)/histbin*1000;
        meanFR    = mean(this_MPH_FR,1);
        semFR     = std(this_MPH_FR,1) / sqrt(size(this_MPH_FR,1));
        patchPdc  = [meanFR fliplr(meanFR) meanFR(1)] + [-semFR fliplr(semFR) -semFR(1)];
        patchX    = [ 1:size(this_MPH_FR,2) size(this_MPH_FR,2):-1:1 1 ] -0.5;
        nPdcPds   = size(this_MPH_FR,1);
        
%         subplot(hsp(1,1)); hold on
% %         bar(linspace(1,size(this_MPH_raster,2),length(histMP)),histMP,...
% %             'FaceColor',colors(irate+1,:),'EdgeColor','none','BarWidth',1)
%         patch(patchX, patchPdc, 'k','EdgeColor','none','FaceAlpha',0.3)
%         plot((1:size(this_MPH_FR,2))-0.5,meanFR,'Color','k','LineWidth',3)    
%         title('Periodic')
%         hold off
        % - - - - - - - - - - - - - - - - 
        
        
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
            
            
            % - - - - - - - - - - - - - - - -
            %         PLOT Trans
            
            histbin = floor(size(this_MPH_raster,2)/30);
            histMP = []; histMP = binspikecounts(mean(this_MPH_raster,1),histbin)/histbin*1000;
            [xb,yb] = stairs(histMP);
            xb = [0; 0; xb; xb(end)+1; xb(end)+1];
            yb = [0; yb(1); yb; yb(end); 0];
            
            meanFR    = mean(this_MPH_FR,1);
            semFR     = std(this_MPH_FR,1) / sqrt(size(this_MPH_FR,1));
            
            subplot(hsp(np-1,1)); hold on
            patch(patchX, patchPdc, 'k','EdgeColor','none','FaceAlpha',0.5)
            patch(patchX, [meanFR fliplr(meanFR) meanFR(1)] + [-semFR fliplr(semFR) -semFR(1)], colors(irate+1,:),'EdgeColor','none','FaceAlpha',0.5)
            plot((1:size(this_MPH_FR,2))-0.5, meanFR, 'Color',colors(irate+1,:),'LineWidth',3)
%             plot(xb*(size(this_MPH_raster,2)/(length(histMP)+1)),yb,...
%                 'Color',colors(ips,:),'LineWidth',2)
            title(sprintf('Following %s (%i pds)',Info.stim_ID_key{ips},size(this_MPH_raster,1)))
            
            % - - - - - - - - - - - - - - - -
            
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
                
                
                % - - - - - - - - - - - - - - - -
                %           PLOT IR
                
                histbin = floor(size(this_MPH_raster,2)/30);
                histMP = []; histMP = binspikecounts(mean(this_MPH_raster,1),histbin)/histbin*1000;
                meanFR    = mean(this_MPH_FR,1);
                semFR     = std(this_MPH_FR,1) / sqrt(size(this_MPH_FR,1));
                
                subplot(hsp(MPH_rate(ii,:).SeqPos,(thisIR==theseIRs)+1)); hold on
                patch(patchX, patchPdc, 'k','EdgeColor','none','FaceAlpha',0.5)
                patch(patchX, [meanFR fliplr(meanFR) meanFR(1)] + [-semFR fliplr(semFR) -semFR(1)], colors(irate+1,:),'EdgeColor','none','FaceAlpha',0.5)
                plot((1:size(this_MPH_FR,2))-0.5, meanFR, 'Color',colors(irate+1,:),'LineWidth',3)
%                 bar(linspace(1,size(this_MPH_raster,2),length(histMP)),histMP,...
%                     'FaceColor',colors(irate+1,:),'EdgeColor','none','BarWidth',1)
                title(sprintf('%s, prevPd: %ih Hz (%i pds)',Info.stim_ID_key{thisIR},ipp,size(this_MPH_raster,1)))
                hold off
                % - - - - - - - - - - - - - - - -
                
            end %PPds
        end %thisIR
        
        
        %==================================================================
        %% Finish and save plot
        
        suptitle(sprintf('%iHz MPH by context, %s type units (%i Pdc pds)',this_rate,UnType{:},nPdcPds))
        
        savedir = fullfile(fn.processed,'MPHplots',SESSION);
        if ~exist(savedir,'dir')
            mkdir([savedir '/svg']);
            mkdir([savedir '/eps']);
        end
        savename = sprintf('MPH_%iHz_%s_%s_all%s',this_rate,SUBJECT,SESSION,UnType{:});
        
%         print_eps_kp(hfe(irate),fullfile([savedir '/eps'],savename))
        
  
    end  %this_rate
    
end %UnType



%% Temporal relationship between unit types
keyboard
%upshot: slightly more common for NS to occur ~1.5 ms before RS

diffs = [];
for ip = 1:length(spiketimes_NS)
    foo = spiketimes_RS-spiketimes_NS(ip);
    diffs = [diffs foo(abs(foo)<100)];
end
figure;
histogram(diffs,-19.5:19.5)




end %function




