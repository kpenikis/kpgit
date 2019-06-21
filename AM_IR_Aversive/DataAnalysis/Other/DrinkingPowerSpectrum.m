function DrinkingPowerSpectrum(optargin)
%
%  pp_plot_rasters(subject, session, [channel, clu] )
%    Plots a raster and psth for each stimulus, separating trials by
%    the preceding stimulus. 
%    If just 2 input variables, will plot all SU and MUs from the session.
%    Uses the TrialData (newer) version of saving stimulus info.
%    Excludes datapoints based on: min Ntrials, min FR. Option to exclude MU.
%
%  KP, 2018-04
%

close all

%!!!!!!!!!!!
SUonly = 1;
%!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!



%% Load Resp data table 

fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'RespStruct_allSU'));
Resp = q.Resp;


%% Prepare figures

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];


hf = figure;
set(hf,'Position',largerect)
hold on



add_y = zeros(1,7);
N_un  = zeros(1,7);

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

        
histbins_rate = logspace(0,1.75,30);
histbins_isi = linspace(2,500,50);




%%

% Remove non-responsive units
[sigUnits,Resp] = identifyResponsiveUnits(Resp);
Resp = Resp(sigUnits);


% Plot one stimulus at a time, without separating previous stimuli

for iUn = 1:numel(Resp)
    
    % Get this unit's info
    subject = Resp(iUn).Subject;
    session = Resp(iUn).Session;
    channel = Resp(iUn).Channel;
    clu     = Resp(iUn).Clu;
    
    % Load data files
    fn = set_paths_directories(subject,session,1);
    if (iUn>1 && ~( strcmp(subject,Resp(iUn-1).Subject) && strcmp(session,Resp(iUn-1).Session) )) || iUn==1
        fprintf('Loading sess %s...\n',session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,Resp(iUn-1).Subject) && strcmp(session,Resp(iUn-1).Session) && channel==Resp(iUn-1).Channel ) )  || iUn==1
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes
    spikes = Spikes.sorted(channel);
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    SpikeStream = zeros(size(SpoutStream));
    SpikeStream(spiketimes) = 1;
    
    SessStart = TrialData.onset(1)+2000;
    SpoutStream(1:SessStart-1) = nan;
    
%     pxON = pwelch(SpikeStream(ONspout),1000);
%     plot(10*log10(pxON),'b')
%     hold on
%     pxOFF = pwelch(SpikeStream(OFFspout),1000);
%     plot(10*log10(pxOFF),'k')
    
    
    ONTOspout  = find(diff(SpoutStream)==1);
    OFFOFspout = find(diff(SpoutStream)==-1);
    
    if SpoutStream(1+find(diff(~isnan(SpoutStream))==1)) == 1
        ONTOspout  = [1+find(diff(~isnan(SpoutStream))==1) ONTOspout];
    else
        OFFOFspout = [1+find(diff(~isnan(SpoutStream))==0) OFFOFspout];
    end
    if SpoutStream(end) == 1
        ONTOspout  = [ONTOspout length(SpoutStream)];
    else
        OFFOFspout = [OFFOFspout length(SpoutStream)];
    end    
    
    [~,startPos] = min([ONTOspout(1) OFFOFspout(1)]);
    
    % If the animal doesn't start on the spout, we have to figure out 
    if startPos~=1, keyboard, end
    
    % For each ON and OFF spout segment
    for ii = 1:numel(ONTOspout)
        
        if ii>numel(ONTOspout) || ii>numel(OFFOFspout)
            continue
        end
        
        % Calculate power spectrum
        if numel(ONTOspout(ii):OFFOFspout(ii))>1000
            pxON(:,ii)  = pwelch(SpikeStream(ONTOspout(ii):OFFOFspout(ii)),1000);
        end
        if ii<numel(ONTOspout) && (numel(OFFOFspout(ii):ONTOspout(ii+1))>1000)
            pxOFF(:,ii) = pwelch(SpikeStream(OFFOFspout(ii):ONTOspout(ii+1)),1000);
        end
        
        % Get Spktime rate distribution
        if numel(ONTOspout(ii):OFFOFspout(ii))>1000
            ONspkRates(ii,:)  = histcounts(1000./diff(find( SpikeStream(ONTOspout(ii):OFFOFspout(ii)) )),histbins_rate);
%             ONspkISIs(ii,:)  = foo.Values;
        end
        if ii<numel(ONTOspout) && (numel(OFFOFspout(ii):ONTOspout(ii+1))>1000)
            OFFspkRates(ii,:) = histcounts(1000./diff(find( SpikeStream(OFFOFspout(ii):ONTOspout(ii+1)) )),histbins_rate);
%             OFFspkISIs(ii,:) = foo.Values;
        end
        
        % Get Spktime rate distribution
        if numel(ONTOspout(ii):OFFOFspout(ii))>1000
            ONspkISIs(ii,:)  = histcounts(diff(find( SpikeStream(ONTOspout(ii):OFFOFspout(ii)) )),histbins_isi);
%             ONspkISIs(ii,:)  = foo.Values;
        end
        if ii<numel(ONTOspout) && (numel(OFFOFspout(ii):ONTOspout(ii+1))>1000)
            OFFspkISIs(ii,:) = histcounts(diff(find( SpikeStream(OFFOFspout(ii):ONTOspout(ii+1)) )),histbins_isi);
%             OFFspkISIs(ii,:) = foo.Values;
        end
        
        
    end
    
    xx = histbins_isi(1:end-1) + mode(diff(histbins_isi))/2;
    
    % Plot average rate distribution from ISIs
    figure(hf); clf;
    patch(1000./[7 7 11 11 7],[0 0.25 0.25 0 0],'y','FaceAlpha',0.2,'EdgeColor','none')
    hold on
    plot(xx,mean(ONspkISIs,1)./sum(mean(ONspkISIs,1)),'b','LineWidth',2)
    plot(xx,mean(OFFspkISIs,1)./sum(mean(OFFspkISIs,1)),'k','LineWidth',2)
    xlim([0 300])
    
%     b=bar([mean(ONspkISIs,1)./sum(mean(ONspkISIs,1));...
%         mean(OFFspkISIs,1)./sum(mean(OFFspkISIs,1))]');
%     b(1).FaceColor = 'c';
%     b(2).FaceColor = 'k';
    
    
%     % Plot average SPECTRUM of each period
%     figure(hf); clf
%     plot(10*log10(mean(pxON,2)),'b','LineWidth',2)
%     hold on
%     xlim([0 50])
%     plot(10*log10(mean(pxOFF,2)),'k','LineWidth',2)
%     
    pause(2)
    
end %iUn



keyboard


keyboard

end