function MPH_PhaseComp(RERUN)
%  MPH_PhaseComp
%   After zoom call with shihab, need to check distribution of mean phases
%   from only steady state AM.
%
% KP, 2020-04
% 


close all


global fn AMrates RateStream trMin

if nargin<1
    RERUN=0;
end

datasource = '250ms'; '1pd'; 


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------
trMin    = 5;

% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',12)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

savedir = fullfile(fn.figs,'PopulationTuning','Phase',datasource);
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%

if RERUN
    
Phs_SS_all       = nan(300,5);
Phs_SS_sig       = nan(300,5);
Phs_First_4      = nan(300,5);
Phs_First_32     = nan(300,5);
Phs_First_4_sig  = nan(300,5);
Phs_First_32_sig = nan(300,5);

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
    
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    
    %% Get MPH data
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes,RateStream);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    %% Calculate mean phase for steady state
    
    for this_rate = AMrates
        
        Period     = 1000/this_rate;
        
        
        switch datasource
            
            %..............................................................
            %..............................................................
            case '1pd'
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Find LAST periodic period (with enough trials)
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~
                iPd        = [];
                ifilt      = MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevPd==this_rate;
                lastPd     = max(MPH.Starttime(ifilt));
                Candidates = MPH(ifilt & MPH.Starttime==lastPd,:);
                [~,iPd]    = min(Candidates.PrevStimID-Candidates.ThisStimID);
                if isempty(iPd)
                    continue
                end
                if numel(iPd)>1
                    keyboard
                end
                
                % Make sure there are at least 10 trials, or else skip
                if size(Candidates(iPd,:).raster{:},1)<10
                    [~,iPd] = max(cellfun(@(x) size(x,1),Candidates.raster));
                end
                if size(Candidates(iPd,:).raster{:},1)<10
                    continue
                end
                
                Spiketimes = Candidates(iPd,:).x{:};
                
                Phs_SS_all(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_SS_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Also get FIRST period, when following 4 and 32 hz
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~
                iPd        = [];
                ifilt      = MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevPd==4;
                iPd        = find(ifilt & MPH.Starttime==0);
                if isempty(iPd)
                    continue
                end
                Spiketimes = MPH(iPd,:).x{:};
                
                Phs_First_4(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_First_4_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
                
                iPd        = [];
                ifilt      = MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevPd==32;
                iPd        = find(ifilt & MPH.Starttime==0);
                if isempty(iPd)
                    continue
                end
                Spiketimes = MPH(iPd,:).x{:};
                
                Phs_First_32(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_First_32_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
            %..............................................................
            %..............................................................
            case '250ms'
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Find LAST 250 ms of periodic stimulus
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                iPds = []; Spiketimes=[];
                
                ifilt    = find(MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevPd==this_rate);
                if this_rate==2
                    iPds = MPH.Starttime(ifilt)==500;
                else
                    iPds = MPH.Starttime(ifilt)>749 & MPH.Starttime(ifilt)<1000;
                end
                
                Data     = MPH(ifilt(iPds),:);
                if isempty(Data) || sum(cellfun(@(x) size(x,1),Data.raster)) < 10
                    continue
                end
                
                Spiketimes = cat(2,Data.x{:});
                
                Phs_SS_all(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_SS_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Also get FIRST period, when following 4 and 32 hz
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % !!! 4 Hz 
                iPds = []; Spiketimes=[];
                
                ifilt    = find(MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevStimID==3);
                if this_rate<=4
                    iPds = MPH.Starttime(ifilt)==0;
                else
                    iPds = MPH.Starttime(ifilt)>=0 & MPH.Starttime(ifilt)<250;
                end
                
                Data     = MPH(ifilt(iPds),:);
                if isempty(Data) || sum(cellfun(@(x) size(x,1),Data.raster)) < 10
                    continue
                end                
                
                Spiketimes = cat(2,Data.x{:});
                
                Phs_First_4(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_First_4_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
                % !!! 32 Hz 
                iPds = []; Spiketimes=[];
                
                ifilt    = find(MPH.AMrate==this_rate & MPH.ThisStimID<7 & MPH.PrevStimID==6);
                if this_rate<=4
                    iPds = MPH.Starttime(ifilt)==0;
                else
                    iPds = MPH.Starttime(ifilt)>=0 & MPH.Starttime(ifilt)<250;
                end
                
                Data     = MPH(ifilt(iPds),:);
                if isempty(Data) || sum(cellfun(@(x) size(x,1),Data.raster)) < 10
                    continue
                end                
                
                Spiketimes = cat(2,Data.x{:});
                
                Phs_First_32(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                
                if UnitData(iUn).VSdata_spk(2,find(this_rate==AMrates)+1)>=13.1
                    Phs_First_32_sig(iUn,this_rate==AMrates) = meanphase(Spiketimes,Period);
                end
                
        end
        
    end %this_rate
    
end %iUn

save(fullfile(savedir,'PhaseData'),'Phs_SS_all','Phs_SS_sig','Phs_First_4','Phs_First_32','Phs_First_4_sig','Phs_First_32_sig','-v7.3')

else
    
    load(fullfile(savedir,'PhaseData'))
    
end


%% POLAR HISTO PLOTS

iRS = UnitInfo.TroughPeak>0.43;
iNS = UnitInfo.TroughPeak<=0.43;
% Fill rest with 0s
iRS(numel(iRS)+1:300) = 0;
iNS(numel(iNS)+1:300) = 0;

degvector = deg2rad(0:22.5:360);

CellType = {'RS' 'NS'};
rmaxvals = [60 24];

for ict = 1:numel(CellType)
    
    figure; hold on
    set(gcf,'Position',fullscreen)
%     suptitle([CellType{ict} ' cells'])
    
    hist_SS_all = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_SS_all);
    hist_SS_sig = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_SS_sig);
    hist_4_all  = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_First_4);
    hist_4_sig  = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_First_4_sig);
    hist_32_all = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_First_32);
    hist_32_sig = makepolarhisto(degvector,eval(['i' CellType{ict}]),Phs_First_32_sig);
    
    
    for ir = 1:numel(AMrates)
        
        subplot(3,5,ir)
        polarplot(degvector,hist_SS_all(:,ir),'Color',[0.7 0.7 1],'LineWidth',3)
        hold on
        polarplot(degvector,hist_SS_sig(:,ir),'Color','b','LineWidth',3)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        
        set(gca,'rlim',[0 rmaxvals(ict)],'rtick',linspace(0,rmaxvals(ict),4),'rticklabel',{'' '' '' num2str(rmaxvals(ict))})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        
        % Stats
        mu = round(rad2deg(meanphase(Phs_SS_sig(eval(['i' CellType{ict}]) & ~isnan(Phs_SS_sig(:,ir)),ir),2*pi)));
        [VS,RS,Pvs] = vectorstrength(Phs_SS_sig(eval(['i' CellType{ict}]) & ~isnan(Phs_SS_sig(:,ir)),ir),2*pi);
        pval = circ_rtest(Phs_SS_sig(eval(['i' CellType{ict}]) & ~isnan(Phs_SS_sig(:,ir)),ir));
        [Pvs; pval];
        
        % kuiper test (bootstrapped)
% %         pval=nan(100,1);
% %         for i=1:100
% %             RndUnDist = round(rand(numel(Phs_SS_sig(~isnan(Phs_SS_sig(:,ir)),ir)),1)*2*pi,4);
% %             pval(i) = circ_kuipertest(round(Phs_SS_sig(~isnan(Phs_SS_sig(:,ir)),ir),4),...
% %                 RndUnDist,200,0);
% %         end
% %         
        
%         title([num2str(AMrates(ir)) ' Hz STEADY STATE'])
        title(sprintf('%i Hz STEADY STATE\np=%0.3f, mu=%i',AMrates(ir),Pvs,mu))
        
        
        subplot(3,5,ir+5)
        polarplot(degvector,hist_4_all(:,ir),'Color',0.7*[1 1 1],'LineWidth',3)
        hold on
        polarplot(degvector,hist_4_sig(:,ir),'Color','k','LineWidth',3)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        
        set(gca,'rlim',[0 rmaxvals(ict)],'rtick',linspace(0,rmaxvals(ict),4),'rticklabel',{'' '' '' num2str(rmaxvals(ict))})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        
        title([num2str(AMrates(ir)) ' Hz, following 4 Hz'])
        
        
        subplot(3,5,ir+10)
        polarplot(degvector,hist_32_all(:,ir),'Color',0.7*[1 1 1],'LineWidth',3)
        hold on
        polarplot(degvector,hist_32_sig(:,ir),'Color','k','LineWidth',3)
        set(gca,'ThetaZeroLocation','bottom',...
            'ThetaDir','clockwise');
        set(gca,'Color','none','RAxisLocation',210)
        
        set(gca,'rlim',[0 rmaxvals(ict)],'rtick',linspace(0,rmaxvals(ict),4),'rticklabel',{'' '' '' num2str(rmaxvals(ict))})
        set(gca,'ThetaTick',0:30:330,'ThetaTickLabel',{'0' '' '' '90' '' ''  '180' '' '' '270' '' ''})
        
        title([num2str(AMrates(ir)) ' Hz, following 32 Hz'])
        
    end %ir
        
    savename = ['MeanPhase_Polar_' CellType{ict}];
    print_eps_kp(gcf,fullfile(savedir,savename))
    
end %celltype




%% 1D HISTO PLOTS

histvec = linspace(0,2*pi,17);


hf_RS=figure; 
set(gcf,'Position',fullscreen)

ymax  = 35;

for ir=1:5
    
    subplot(3,5,ir)
    histogram(Phs_SS_all(iRS,ir),histvec,'FaceColor','b','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_SS_sig(iRS,ir),histvec,'FaceColor','b','FaceAlpha',1,'EdgeColor','none')
%     xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz STEADY STATE'])
    
    
    subplot(3,5,ir+5)
    histogram(Phs_First_4(iRS,ir),histvec,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_First_4_sig(iRS,ir),histvec,'FaceColor','k','FaceAlpha',1,'EdgeColor','none')
%     xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz, following 4 Hz'])
    
    
    subplot(3,5,ir+10)
    histogram(Phs_First_32(iRS,ir),histvec,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_First_32_sig(iRS,ir),histvec,'FaceColor','k','FaceAlpha',1,'EdgeColor','none')
    xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz, following 32 Hz'])
    
end
suptitle('RS cells')

savename = 'MeanPhaseDistr_RScells';
print_eps_kp(gcf,fullfile(savedir,savename))



hf_NS=figure; 
set(gcf,'Position',fullscreen)

ymax  = 13;

for ir=1:5
    
    subplot(3,5,ir)
    histogram(Phs_SS_all(iNS,ir),histvec,'FaceColor','b','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_SS_sig(iNS,ir),histvec,'FaceColor','b','FaceAlpha',1,'EdgeColor','none')
%     xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz STEADY STATE'])
    
    subplot(3,5,ir+5)
    histogram(Phs_First_4(iNS,ir),histvec,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_First_4_sig(iNS,ir),histvec,'FaceColor','k','FaceAlpha',1,'EdgeColor','none')
%     xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz, following 4 Hz'])
    
    subplot(3,5,ir+10)
    histogram(Phs_First_32(iNS,ir),histvec,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','none')
    hold on
    histogram(Phs_First_32_sig(iNS,ir),histvec,'FaceColor','k','FaceAlpha',1,'EdgeColor','none')
    xlabel('Mean phase (radians)')
    ylabel('N cells')
    xlim([0 2*pi])
    ylim([0 ymax])
    title([num2str(AMrates(ir)) ' Hz, following 32 Hz'])
    
end
suptitle('NS cells')

savename = 'MeanPhaseDistr_NScells';
print_eps_kp(gcf,fullfile(savedir,savename))


keyboard


end %function




