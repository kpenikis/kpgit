function ap_IRtrans_allRasters
%
%  ap_IRtrans_diffspks_v2(select_subject,select_session,channels,clus)
%
%  KP, 2018-01
 
close all


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 2*scrsz(4)/3 scrsz(3) 2*scrsz(4)/3];

histbinsize = 5;
smthbinsize = 50;


colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors = [colors; 0.7.*bone(4)];

pbcolors = colors(3:4,:);
% pbcolors  = [  0 100   0;...
%               16  78 139]./255;

% wincolors = flipud(copper(numel(t_win)-1));
wincolors = flipud(winter(numel(t_win)-1));



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fn = set_paths_directories;
T = readtable([fn.processed '/Transitions/IRtransitionsCDF_SU.txt']);
Tsig250 = sortrows(T(T.zScore_250>2,:),'zScore_250','descend');


for iUn = 1:size(Tsig250,1)
    
    % Get datapoint info
    subject = Tsig250.Subject{iUn};
    session = Tsig250.Session{iUn};
    channel = Tsig250.ch(iUn);
    clu     = Tsig250.clu(iUn);
    spl     = Tsig250.stimpars_1(iUn);
    lpn     = Tsig250.stimpars_2(iUn);
    amd     = 1;
    ib      = Tsig250.IRblock(iUn);
    
    % Load necessary files
    fn = set_paths_directories(subject,session);
    fprintf('Session %s:  loading data...\n',session)
    filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_SoundData',subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    % Get spiketimes, artifact, and fs
    spikes = Spikes.sorted(channel);
    spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
    ArtifactFlag = Info.artifact(channel).SDsamples;
    fs = round(Info.fs_sound);
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    
    % Get smoothed FR / convert to z-score
    
    Stream_Spks = zeros(1,1000*ceil((size(SoundData,2)/Info.fs_sound)));
    Stream_Spks(spiketimes) = 1;
    
    %either with standard 20 ms bin
    Stream_FRbin = 1000*(binspikecounts(Stream_Spks,histbinsize)/histbinsize);
    Stream_FRbin(isinf(Stream_FRbin)) = nan;
    foo = repmat(Stream_FRbin,histbinsize,1);
    Stream_FR = reshape(foo,1,histbinsize*length(Stream_FRbin));
    Stream_FR = Stream_FR(1:ceil(length(SoundData)/Info.fs_sound*1000));
    
    %or with sliding X ms boxcar
    Stream_FRsmooth = smoothFR(Stream_Spks,smthbinsize);
    Stream_FRsmooth = Stream_FRsmooth(1:ceil(length(SoundData)/Info.fs_sound*1000));
    
    % Determine the time of the beginning of actual recording session
    % max of first spiketime or 5 seconds before unmodulated sound came
    % on
    sampStart = find(diff(SoundData(8,:))==-1);
    msStart   = max( spiketimes(1), round(sampStart(1)/Info.fs_sound*1000)-5000 );
    
    % Convert FR to z-score
    Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
    Stream_zscore = [nan(1,msStart-1) Stream_zscore];
    
    % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    
    fprintf('plotting ch %i clu %i...\n',channel,clu)
    
    % Get this unit/stimulus combo's N clean blocks (trials)
    [~,minTime] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd,7:10);
    t_win = 1:tstep:minTime;
    t1 = -2199;
    
    % Make figure
    hfex = figure;
    set(hfex,'Position',figsize1,'NextPlot','add')
    hold on
    
    % Get this block start times
    [bkStart_samps,~,...
        bkStart_ms,~]  =  get_blockOnsets( SoundData,...
        ib,spl,lpn,amd,ArtifactFlag,fs);
    
    % Get list of blocks that came before this one
    prevBlocks = sort(unique(SoundData(8,bkStart_samps-1)));
    
    % Set up some empty vars
    legstr = cell(1,numel(prevBlocks)); clear ip
    raster_x_1 = [];
    raster_y_1 = [];
    raster_x_2 = [];
    raster_y_2 = [];
    stim = nan( numel(prevBlocks), minTime-t1+1, ceil(numel(bkStart_ms)/2) );
    psth = nan( numel(prevBlocks), minTime-t1+1, ceil(numel(bkStart_ms)/2) );
    IR_psth_wins = nan( numel(prevBlocks), mode(diff(t_win)), numel(t_win)-1 );
    
    for pb = prevBlocks
        
        ii = find(SoundData(8,bkStart_samps-1) == pb);
        
        % Collect FR for specified transitions
        for jj = 1:numel(ii)
            psth(pb==prevBlocks,:,jj) = Stream_FRsmooth( (bkStart_ms(ii(jj))+t1) : ((bkStart_ms(ii(jj))+minTime)) );
            stim(pb==prevBlocks,:,jj) = envelope(SoundData(2, round(( (bkStart_ms(ii(jj))+t1) : ((bkStart_ms(ii(jj))+minTime)))  ./1000.*Info.fs_sound) ),20,'rms');
            sp=[]; sp = spiketimes( spiketimes>=(bkStart_ms(ii(jj))+t1) & spiketimes<((bkStart_ms(ii(jj))+minTime)) ) - bkStart_ms(ii(jj))+1;
            eval( sprintf('raster_x_%i = [raster_x_%i sp];',find(pb==prevBlocks),find(pb==prevBlocks)) )
            eval( sprintf('raster_y_%i = [raster_y_%i jj*ones(1,numel(sp))];',find(pb==prevBlocks),find(pb==prevBlocks)) )
        end
        
        % Save the psth data for each time window
        for iw = 2:numel(t_win)
            IR_psth_wins(pb==prevBlocks,:,iw-1) = mean(psth(pb==prevBlocks,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
        end
        
        
    end %pb
    
    
    for iw = 2:numel(t_win)
        % Plot the space between PSTHs
        hs(3)=subplot(5,1,4:5); hold on
        fill( [t_win(iw-1)+1:t_win(iw) t_win(iw):-1:t_win(iw-1)+1] , [IR_psth_wins(1,:,iw-1) fliplr(IR_psth_wins(2,:,iw-1))] , wincolors(iw-1,:),...
            'EdgeColor','none','FaceAlpha',0.65);
    end
    
    % Now plot each PSTH  and raster and stim rms trace
    add_y = 0;
    for pb = prevBlocks
        subplot(hs(3)); hold on; box off
        ip(pb==prevBlocks) = plot(t1:minTime,mean(psth(pb==prevBlocks,:,:),3,'omitnan'),'Color',pbcolors(pb==prevBlocks,:),'LineWidth',4);
        hold on
        legstr{pb==prevBlocks} = [Info.blockKey{pb} ' Hz, n=' num2str(jj)];
        
        hs(1)=subplot(5,1,1); hold on; box off
        plot(t1:minTime,mean(stim(pb==prevBlocks,:,:),3,'omitnan'),'Color',pbcolors(pb==prevBlocks,:),'LineWidth',4)
        set(gca,'xtick',[],'ytick',[])
        
        hs(2)=subplot(5,1,2:3); hold on; box off
        plot(eval(['raster_x_' num2str(find(pb==prevBlocks))]),eval(['raster_y_' num2str(find(pb==prevBlocks))])+add_y,...
            '.','MarkerSize',15,'Color',pbcolors(pb==prevBlocks,:))
        add_y = add_y + max(eval(['raster_y_' num2str(find(pb==prevBlocks))]));
        
    end
    
    % Finish plot settings
    linkaxes(hs,'x')
    xlim([-2100 minTime])
    
    subplot(hs(2))
    set(gca,'ylim',[0 add_y+1],'xtick',[])
    
    subplot(hs(3))
    ylim([0 10*ceil(max(max(max(IR_psth_wins)))/10)+5])
    xlabel('Time from transition (ms)')
    ylabel('Spikes/s')
    legend(ip,legstr)
    
    suptitle(sprintf('Transition pdc to IR_%i   |   %s %s ch%i clu%i   |   %idB SPL, 100-%i Hz   |   zscore 0-250 ms = %0.3f',ib,subject,session,channel,clu,spl,lpn,Tsig250.zScore_250(iUn)))
    
    % Save this plot
    savedir = fullfile(fn.processed,'Transitions','SigUnitRasters');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(hfex,fullfile(savedir,sprintf('IRtransition_ranked%i_%s_%i_%i_bk%i',iUn,session,channel,clu,ib)))
    
    
    aaa=235;
    
end %iUn



end %function




