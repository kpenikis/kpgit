function OnsetWinRatio( UnitData, UnitInfo )
% Run after AssessUnits. Make sure to manually add unit strings to
% SplitUnits mat file first. Skips units that have already been merged.
% 
%  KP, 2018-06
%



close all
global fn AMrates rateVec_AC rateVec_DB

degrees = [15 75];
winBeg = deg2rad(degrees(1));
winEnd = deg2rad(degrees(2));

rng('shuffle')


%% Load Unit data files

fn = set_paths_directories('','',1);

if nargin<1
    q = load(fullfile(fn.processed,'Units'));
    UnitData = q.UnitData;
    UnitInfo = q.UnitInfo;
    clear q
end

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];



%% Settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];


% % Prepare figure
% hf = figure;
% set(hf,'Position',fullscreen)


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
        


%% 

% thisUnit = find(strcmp(UnitInfo.Session,'QB') & UnitInfo.Channel==3 & UnitInfo.Clu==32);


OWRatio_BMF = nan(numel(UnitData),1);
OWRatio_BVS = nan(numel(UnitData),1);
OWR_2hz     = nan(numel(UnitData),1);

for iUn = 1:numel(UnitData)
    
    %%% still must edit for merged units
    if strncmp(UnitInfo.RespType{iUn},'merged',6)
        
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
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
        
        % Get spiketimes and shift based on calculated integration time
        spiketimes = round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000);  %ms
        spiketimes = spiketimes-round(UnitData(iUn).IntTime_spk);
        
        % Get sound parameters
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        MPH2 = makeMPHtable(TrialData,Info.artifact(channel).trials,dBSPL,LP,spiketimes);
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        % Merge the tables
        MPH = [MPH; MPH2];
        
    end
    
    
    
    
    %% Get data for each AM rate
    
    
    hf = figure;
    set(hf,'Position',fullscreen)
    plot([1 64],[0 0],'--k')
    hold on
    set(gca,'xscale','log','xlim',[1 64],'xtick',AMrates)
    
    OWRatio = nan(3,numel(AMrates));
    
    for irate = AMrates
        
        Theta       =  linspace(0,2*pi,ceil(1000/irate));
        gaussWinLen =  round(30/360*1000/irate); %30 degrees converted to ms
        
        idx      =  find(MPH.AMrate==irate);
        idx_pdc  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID<7  & MPH(MPH.AMrate==irate,:).PrevPd==irate)';
        idx_ira  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==7 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        idx_irb  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==8 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        
        FinalData = struct;
        
        
        % Periodic context instances of this rate
        
        excl_idx = MPH(idx_pdc,:).PrevStimID==(1+find(irate==AMrates)) & MPH(idx_pdc,:).Starttime<100;
        idx_pdc(excl_idx) = [];
        FinalData(1).PrevPd = irate;
        FinalData(1).raster = cell2mat(MPH(idx_pdc,:).raster);
        
        raster   =  mean(cell2mat(MPH(idx_pdc,:).raster),1);
        OWRatio(1,irate==AMrates)  =  ( mean(raster(:,(Theta>=winBeg)&(Theta<=winEnd))) - mean(raster) ) / mean(raster);
        
        
        % Irregular sequence A
        
        % Combine rows that have the same previous period
        unqPrPds = unique(MPH(idx_ira,:).PrevPd)';
        for ipp = unqPrPds
            theseidx = idx_ira(MPH(idx_ira,:).PrevPd==ipp);
            FinalData(end+1).PrevPd = ipp;
            FinalData(end).raster   = cell2mat( MPH(theseidx,:).raster );
        end
        
        raster   =  mean(cell2mat(MPH(idx_ira,:).raster),1);
        if ~isempty(raster)
            OWRatio(2,irate==AMrates)  =  ( mean(raster(:,(Theta>=winBeg)&(Theta<=winEnd))) - mean(raster) ) / mean(raster);
        end
        
        
        % Irregular sequence B
        
        % Combine rows that have the same previous period
        unqPrPds = unique(MPH(idx_irb,:).PrevPd)';
        for ipp = unqPrPds
            theseidx = idx_irb(MPH(idx_irb,:).PrevPd==ipp);
            FinalData(end+1).PrevPd = ipp;
            FinalData(end).raster   = cell2mat( MPH(theseidx,:).raster );
        end
        
        raster   =  mean(cell2mat(MPH(idx_irb,:).raster),1);
        if ~isempty(raster)
            OWRatio(3,irate==AMrates)  =  ( mean(raster(:,(Theta>=winBeg)&(Theta<=winEnd))) - mean(raster) ) / mean(raster);
        end
        
        
        
        %% Calculate ratio for each trial, with separated prev pds
        
        minTrs = min(cellfun(@(x) size(x,1),{FinalData.raster}));
        for ic = 1:numel(FinalData)
            rndTrs = randperm(size(FinalData(ic).raster,1));
            data = FinalData(ic).raster(rndTrs(1:minTrs),:);
            FinalData(ic).OWR = ( mean(data(:,(Theta>=winBeg)&(Theta<=winEnd)),2) - mean(data,2) ) ./ mean(data,2);
            plot( irate-irate*0.05+irate*0.1*rand(size(FinalData(ic).OWR)) ,...
                FinalData(ic).OWR,'.','Color',colors(1+find(FinalData(ic).PrevPd==AMrates),:))
        end
        
        
        
        
    end %irate
    
    % Save Onset Window Ratio for BMF stimuli
    UnitData(iUn).OWR = OWRatio;
    if ~isempty(UnitData(iUn).iBMF_FR)
        OWRatio_BMF(iUn) = OWRatio(1,UnitData(iUn).iBMF_FR) ;
    end
    if ~isempty(UnitData(iUn).iBMF_VS)
        OWRatio_BVS(iUn) = OWRatio(1,UnitData(iUn).iBMF_VS) ;
    end
    
    
    %% Plot as tuning curve
    
%     hf = figure; 
%     plot([1 64],[0 0],'--k')
%     hold on
    plot(AMrates,OWRatio(1,:),'o-k','MarkerSize',15,'LineWidth',3)
    plot(AMrates,OWRatio(2,:),'o-','Color',colors(7,:),...
        'MarkerSize',10,'LineWidth',2)
    plot(AMrates,OWRatio(3,:),'o-','Color',colors(8,:),...
        'MarkerSize',10,'LineWidth',2)
%     set(gca,'xscale','log','xlim',[1 64],'xtick',AMrates)
    xlabel('AM rate of cycle')
    ylabel('Onset Window Ratio')
    title(sprintf('%s %s ch%i clu%i  |  RespType: %s\n%i-%i degrees',subject,session,channel,clu,UnitInfo.RespType{iUn},degrees(1),degrees(2)))
    
    
    %% Save figure
    
    savedir = fullfile(fn.processed,'OWRatio','Units');
    if ~exist(savedir)
        mkdir(fullfile(savedir,'eps'))
        mkdir(fullfile(savedir,'svg'))
    end
    
    savename = sprintf('OWRtuning_%s_%s_%i_%i_win%ito%i',subject,session,channel,clu,degrees(1),degrees(2));
    print_eps_kp(hf,fullfile(savedir,'eps',savename))
    print_svg_kp(hf,fullfile(savedir,'svg',savename))
    
    close(hf)
    
    
    %% Save 2 Hz OWR to store in UnitInfo table
    
    OWR_2hz(iUn,1) = OWRatio(1,1);
    
    
end %iUn


% Store 2 Hz OWR in UnitInfo table
UnitInfo.OWR_2hz = OWR_2hz;
keyboard %to re-save, if desired


%% Plot histograms of OWR values across all units

OWRs    = vertcat(UnitData.OWR);
BaseFRs = vertcat(UnitData.BaseFR);

histvec = [-1:0.25:5];

figure;
subplot(2,3,1)
histogram(OWRs(:,1),histvec)
subplot(2,3,2)
histogram(OWRs(:,2),histvec)
subplot(2,3,3)
histogram(OWRs(:,3),histvec)
subplot(2,3,4)
histogram(OWRs(:,4),histvec)
subplot(2,3,5)
histogram(OWRs(:,5),histvec)

subplot(2,3,6)
[~,imx] = max(abs(OWRs),[],2);
histogram(OWRs(sub2ind(size(OWRs),1:size(OWRs,1),imx')),histvec)



figure;

subplot(1,3,1)
histogram(OWRatio_BMF,histvec)

subplot(1,3,2)
histogram(OWRatio_BVS,histvec)

subplot(1,3,3)
[~,imx] = max(abs(OWRs),[],2);
histogram(OWRs(sub2ind(size(OWRs),1:size(OWRs,1),imx')),histvec)





end
