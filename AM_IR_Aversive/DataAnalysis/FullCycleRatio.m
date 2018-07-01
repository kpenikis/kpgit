function FullCycleRatio( UnitData, UnitInfo )
% Run after AssessUnits. Make sure to manually add unit strings to
% SplitUnits mat file first. Skips units that have already been merged.
% 
%  KP, 2018-06
%



close all
global fn AMrates rateVec_AC rateVec_DB

winWidth = 45;
winStep  = 15;
Wins = [ (0:winStep:359)' winWidth+(0:winStep:359)'];
Theta = deg2rad(mean(Wins,2));

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


for ir = 1:numel(AMrates)
    hf(ir) = figure;
    set(hf(ir),'Position',fullscreen)
    hs(ir) = subplot(1,1,1,polaraxes);
    polarplot(linspace(0,2*pi,100),ones(1,100),'g-')
    ax=gca;
    ax.ThetaZeroLocation = 'bottom';
    ax.ThetaDir = 'clockwise';
    rlim([0 6])
    title([num2str(AMrates(ir)) ' Hz'])
    
    hfc(ir) = figure;
    set(hfc(ir),'Position',fullscreen)
    hsc(ir) = subplot(1,1,1,polaraxes);
    polarplot(linspace(0,2*pi,100),ones(1,100),'g-')
    ax=gca;
    ax.ThetaZeroLocation = 'bottom';
    ax.ThetaDir = 'clockwise';
    rlim([0 6])
    title([num2str(AMrates(ir)) ' Hz'])
end


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
        
colSparse    = [0.157 0.44  1]; %[0.208 0.678 0.275];
colSustained = [1     0.44  0.157]; %[0.200 0.282 0.800];
colGap       = [0.139 0.159 0.190];


%% 

% thisUnit = find(strcmp(UnitInfo.Session,'QB') & UnitInfo.Channel==3 & UnitInfo.Clu==32);


OWRatio_BMF = nan(numel(UnitData),1);
OWRatio_BVS = nan(numel(UnitData),1);


for iUn = 1:numel(UnitData)
    
    %%
    if strncmp(UnitInfo.RespType{iUn},'merged',6)
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
    
    
    OWRatio = nan(size(Wins,1),numel(AMrates),3);
    
    for irate = AMrates
        
        TimeRadians  =  linspace(0,2*pi,ceil(1000/irate));
        
        idx      =  find(MPH.AMrate==irate);
        idx_pdc  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID<7  & MPH(MPH.AMrate==irate,:).PrevPd==irate)';
        idx_ira  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==7 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        idx_irb  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==8 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        
        
        for iw = 1:size(Wins,1)
            
            winBeg   = deg2rad(mod(Wins(iw,1),360));
            winEnd   = deg2rad(mod(Wins(iw,2),360));
            
            if winEnd>winBeg
                Window = (TimeRadians>=winBeg)&(TimeRadians<=winEnd);
            else
                Window = (TimeRadians>=winBeg)|(TimeRadians<=winEnd);
            end
            
            % Periodic context instances of this rate
            excl_idx = MPH(idx_pdc,:).PrevStimID==(1+find(irate==AMrates)) & MPH(idx_pdc,:).Starttime<100;
            idx_pdc(excl_idx) = [];
            
            raster   =  mean(cell2mat(MPH(idx_pdc,:).raster),1);
            OWRatio(iw,irate==AMrates,1)  =  ( mean(raster(:,Window)) - mean(raster) ) / mean(raster);
            
            
            % Irregular sequence A
            raster   =  mean(cell2mat(MPH(idx_ira,:).raster),1);
            if ~isempty(raster)
                OWRatio(iw,irate==AMrates,2)  =  ( mean(raster(:,Window)) - mean(raster) ) / mean(raster);
            end
            
            
            % Irregular sequence B
            raster   =  mean(cell2mat(MPH(idx_irb,:).raster),1);
            if ~isempty(raster)
                OWRatio(iw,irate==AMrates,3)  =  ( mean(raster(:,Window)) - mean(raster) ) / mean(raster);
            end
            
            
        end %iwin
        
        
        %% Plot datapoints
        
        % First try categorizing units by ratios at certain phases
        if (OWRatio(2,irate==AMrates,1)+1)<1 && mean(OWRatio(10:13,irate==AMrates,1)+1)>1
            Label = 'gap';
            thiscolor = colGap;
        elseif (OWRatio(2,irate==AMrates,1)+1)>1 && (OWRatio(6,irate==AMrates,1)+1)<1 && mean(OWRatio(10:13,irate==AMrates,1)+1)<1
            Label = 'sparse';
            thiscolor = colSparse;
        elseif (OWRatio(2,irate==AMrates,1)+1)>1 && (OWRatio(6,irate==AMrates,1)+1)>1 && (OWRatio(6,irate==AMrates,1)+1)<(OWRatio(2,irate==AMrates,1)+1) %&& mean(OWRatio(10:13,irate==AMrates,1)+1)>1
            Label = 'sustained';
            thiscolor = colSustained;
        else
%             keyboard
            thiscolor = 0.7*[1 1 1];
        end
        
        figure(hfc(irate==AMrates));
        subplot(hsc(irate==AMrates)); hold on
        polarplot([Theta; Theta(1)],[OWRatio(:,irate==AMrates,1); OWRatio(1,irate==AMrates,1)]+1,'-','Color',thiscolor)
        
        
        % Then plot example units according to manual caegorization
        
        switch UnitInfo.RespType{iUn}
            case 'sparse'
                thiscolor = colSparse;
            case 'sustained'
                thiscolor = colSustained;
            case 'gap'
                thiscolor = colGap;
            otherwise
                continue
                thiscolor = 0*[1 1 1];
        end
        
        figure(hf(irate==AMrates));
        subplot(hs(irate==AMrates)); hold on
        polarplot([Theta; Theta(1)],[OWRatio(:,irate==AMrates,1); OWRatio(1,irate==AMrates,1)]+1,'-','Color',thiscolor)
        
        
    end %irate
    
    
end %iUn



%% Save figures

savedir = fullfile(fn.processed,'OWRatio','FullCycles');
if ~exist(savedir)
    mkdir(fullfile(savedir,'eps'))
    mkdir(fullfile(savedir,'svg'))
end

for ir = 1:numel(AMrates)
    savename = sprintf('WinRatios_exUnits_%ihz',AMrates(ir));
    print_eps_kp(hf(ir),fullfile(savedir,'eps',savename))
    print_svg_kp(hf(ir),fullfile(savedir,'svg',savename))
    
    savename = sprintf('WinRatios_allUnits_%ihz',AMrates(ir));
    print_eps_kp(hfc(ir),fullfile(savedir,'eps',savename))
    print_svg_kp(hfc(ir),fullfile(savedir,'svg',savename))
end


return




%% Save Unit files again
save(fullfile(fn.processed,'Units'),'UnitInfo','UnitData','-v7.3');


end
