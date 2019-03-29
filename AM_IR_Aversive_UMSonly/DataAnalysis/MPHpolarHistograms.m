function MPHpolarHistograms



close all
global fn AMrates rateVec_AC rateVec_DB

% gaussWinLen = 10;


%% Settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
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



%% 

% thisUnit = find(strcmp(UnitInfo.Session,'QB') & UnitInfo.Channel==3 & UnitInfo.Clu==32);

for iUn = 1:numel(UnitData)
    
    %%% still must edit for merged units
    if strncmp(UnitInfo.RespType{iUn},'merged',6)
        keyboard
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
    
    
    
    %% Prepare figure 
    
    hf = figure;
    set(hf,'Position',fullscreen)
    
    % Set ymaxval
    ceil(2.5*max(MPH.Prev100msFR));
    ymaxval = ceil(mean(MPH.Prev500msFR)+2*range(MPH.Prev100msFR));
    
    % Get data for each period and plot
    for irate = AMrates
        
        Theta       =  linspace(0,2*pi,ceil(1000/irate));
        gaussWinLen =  round(30/360*1000/irate); %30 degrees converted to ms
        
        idx      =  find(MPH.AMrate==irate);
        idx_pdc  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID<7  & MPH(MPH.AMrate==irate,:).PrevPd==irate)';
        idx_ira  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==7 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        idx_irb  =  idx(MPH(MPH.AMrate==irate,:).ThisStimID==8 & MPH(MPH.AMrate==irate,:).PrevPd~=irate)';
        
        % Periodic context instances of this rate
        MPHhistos = cell2mat(cellfun(@mean, MPH(idx_pdc,:).raster,'UniformOutput',false)).*1000;
        MPHhistos = convolveGauss(MPHhistos,gaussWinLen);
        FinalData = combineSamePrevPdHistos(MPHhistos,MPH(idx_pdc,:));
        
        subplot(3,5,5*0+find(irate==AMrates),polaraxes)
        hold on
        for ipp = AMrates
            polarplot([Theta 0],[FinalData(ipp==AMrates,:) FinalData(ipp==AMrates,1)],'LineWidth',3,'Color',colors(find(ipp==AMrates)+1,:))
        end
        rlim([0 ymaxval])
        ax=gca;
        ax.ThetaZeroLocation = 'bottom';
        ax.ThetaDir = 'clockwise';
        title([num2str(irate) ' Hz, pdc'])
                
        % Irregular context "A" instances of this rate
        MPHhistos  = cell2mat(cellfun(@mean, MPH(idx_ira,:).raster,'UniformOutput',false)).*1000;
        if ~isempty(MPHhistos)
            MPHhistos = convolveGauss(MPHhistos,gaussWinLen);
            FinalData = combineSamePrevPdHistos(MPHhistos,MPH(idx_ira,:));
            
            subplot(3,5,5*1+find(irate==AMrates),polaraxes)
            hold on
            for ipp = AMrates
                polarplot(Theta,FinalData(ipp==AMrates,:),'LineWidth',3,'Color',colors(find(ipp==AMrates)+1,:))
            end
            rlim([0 ymaxval])
            ax=gca;
            ax.ThetaZeroLocation = 'bottom';
            ax.ThetaDir = 'clockwise';
            title([num2str(irate) ' Hz, IRa'])
        end
        
        % Irregular context "B" instances of this rate
        MPHhistos  = cell2mat(cellfun(@mean, MPH(idx_irb,:).raster,'UniformOutput',false)).*1000;
        if ~isempty(MPHhistos)
            MPHhistos = convolveGauss(MPHhistos,gaussWinLen);
            FinalData = combineSamePrevPdHistos(MPHhistos,MPH(idx_irb,:));
            
            subplot(3,5,5*2+find(irate==AMrates),polaraxes)
            hold on
            for ipp = AMrates
                polarplot(Theta,FinalData(ipp==AMrates,:),'LineWidth',3,'Color',colors(find(ipp==AMrates)+1,:))
            end
            rlim([0 ymaxval])
            ax=gca;
            ax.ThetaZeroLocation = 'bottom';
            ax.ThetaDir = 'clockwise';
            title([num2str(irate) ' Hz, IRb'])
        end
        
    end %irate
    
    suptitle(sprintf('%s %s %i %i, RespType: %s\n\n',subject,UnitInfo(iUn,:).Session{:},channel,clu,UnitInfo.RespType{iUn}))
    
    
    
    %% Save figure
    
    savedir = fullfile(fn.processed,'PolarPlots',UnitInfo.RespType{iUn});
    if ~exist(savedir,'dir')
        mkdir([savedir '/eps'])
        mkdir([savedir '/svg'])
    end
    savename = sprintf('PolarPlots_%s_%s_%i_%i',subject,UnitInfo(iUn,:).Session{:},channel,clu);
    
    print_eps_kp(hf,fullfile(savedir,'eps',savename))
    print_svg_kp(hf,fullfile(savedir,'svg',savename))
    
    close(hf)
    
    
    
end %iUn

aaa=2342;  



end