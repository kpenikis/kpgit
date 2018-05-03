function ap_IRtrans_diffspks_v2(select_subject,select_session,channels,clus)
%
%  ap_zscore_plots(subject, session, channel, clu)
%    Similar to ap_plot_rasters, but converts FR to a zscore.
%
%  KP, 2017-08
 
% close all

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  20; 
%!!!!!!!!!!!!!!!!!
tstep = 250;
t_win = 1:tstep:1970;
%!!!!!!!!!!!!!!!!!
PLOT_EX = 1;
%!!!!!!!!!!!!!!!!!
if PLOT_EX
    numiterations = 10;
else
    numiterations = 100000;
end
%!!!!!!!!!!!!!!!!!
PLOT_Distr = 1;
%!!!!!!!!!!!!!!!!!


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

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


unType = {'' 'SU' 'MU'};

N = 0;



    
Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.IRblock  = 0;
Zdata.stimpars = [0 0];



%%  SUBJECTS

if nargin<1
    subjects = {'WWWr_244300' 'WWWf_244303'};
else
    subjects = {select_subject};
end

for subj = 1:numel(subjects)
    
    subject = subjects{subj};
    
    
    switch subject
        case 'WWWf_244303'
            region = 'caudal A1';
            subjcol = 'k';
        case 'WWWr_244300'
            region = 'DP/VP';
            subjcol = 'b';
    end
    

%%  SESSIONS

if ~exist('select_session','var')
    
    % Get list of sessions to check for sorted data
    
    fn = set_paths_directories(subject);
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
    % Get list of all sessions with Spikes file
    Sessions = cell(1,1);
    for ifn = 1:numel(SpkFns)
        splitStr = regexp(SpkFns(ifn).name,'_','split');
        splitStr2 = regexp(splitStr{3},'-','split');
        if length(splitStr2{2})==2
            Sessions{end+1,1} = splitStr2{2};
        end
    end
    Sessions = Sessions(2:end,:);
    
else
    Sessions = {select_session};
end


% Step through each session
for sess = Sessions'
    
session = sess{:};

%%
% Load data files
fn = set_paths_directories(subject,session);
fprintf('Session %s:  loading data...\n',session)
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_SoundData',subject,session); load(fullfile(fn.processed,subject,filename));
load(fullfile(fn.stim,'IRsequences.mat'))


%%
% GET STIM INFO

% Get unique dBSPLs
dBSPL = unique(SoundData(4,:));
rm_i=[];
for ii = 1:numel(dBSPL)
    if (numel(find(SoundData(4,:)==dBSPL(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
dBSPL(rm_i) = [];

% Get unique noisebands (based on LP)
LP = unique(SoundData(6,:));
rm_i=[];
for ii = 1:numel(LP)
    if (numel(find(SoundData(6,:)==LP(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
LP(rm_i) = [];

% Get unique AM depths
AMdepth = unique(SoundData(3,:));
rm_i=[];
for ii = 1:numel(AMdepth)
    if (numel(find(SoundData(3,:)==AMdepth(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
AMdepth(rm_i) = [];



%% GET SPIKE TIMES

% Step through all channels if not specified
if nargin<3 && ~exist('channels','var')
    channels = [1:7 9:16];
end


%% STEP THROUGH EACH CHANNEL
for channel = channels
        
    % Artifact for this channel
    ArtifactFlag = Info.artifact(channel).SDsamples;
    fs = round(Info.fs_sound);
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if all(spikes.labels(:,2)==1)
            continue
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = clus'
                
        % !! Only SU for now !!
        if SUonly && (spikes.labels(spikes.labels(:,1)==clu,2) ~= 2)
            continue
        end
        
        try
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        % spiketrials = spikes.trials(unit_in);
        
        if isempty(spiketimes)
            error('no spike events found for this clu')
            %also skip clus with few events 
        elseif numel(spiketimes) < round(FRcutoff*length(SoundData)/Info.fs_sound)
            continue
        end
        
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
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
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        % Convert FR to z-score
        Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
        Stream_zscore = [nan(1,msStart-1) Stream_zscore];
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        
        %%        
        
        % Also step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                for amd = AMdepth
                    
                    close all 
                    %%
                    
                    fprintf('plotting ch %i clu %i...\n',channel,clu)

                    % Get this unit/stimulus combo's N clean blocks (trials)
                    [blocksN,minTime] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd,7:10);
                    
                    t_win = 1:tstep:minTime;
                    
                    %~~~~~~~~~~~
                    t1 = -999;
                    %~~~~~~~~~~~
                                        
                    
                    % Go through each stiulus
                    for ib = 7:10
                        
                        if PLOT_EX
                            RNDsaved = 0;
                            % Make figure
                            hfex = figure;
                            set(hfex,'Position',figsize1,'NextPlot','add')
                            hold on
                        end
                        
                        if blocksN(ib)<minTrs, continue, end
                        
                        % Get this block start times
                        [bkStart_samps,bkStop_samps,...
                            bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
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
                            
                            for iw = 2:numel(t_win)
                                % Save the psth data for each time window
                                IR_psth_wins(pb==prevBlocks,:,iw-1) = mean(psth(pb==prevBlocks,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
                            end %iw
                            
                            
                        end %pb
                        
                        
                        
                        % Calculate the difference between the PSTHs for
                        % each time window
                        TestDiffSpikeCount = nan(1,numel(t_win)-1);
                        for iw = 2:numel(t_win)
                            % Get cumulative count of difference in n
                            % spikes
                            TestDiffSpikeCount(1,iw-1) = sum(abs( diff(IR_psth_wins(:,:,iw-1)./1000,1) ));
                            
                            if PLOT_EX
                                hs(3)=subplot(5,1,4:5); hold on
                                % Plot the space between PSTHs
                                fill( [t_win(iw-1)+1:t_win(iw) t_win(iw):-1:t_win(iw-1)+1] , [IR_psth_wins(1,:,iw-1) fliplr(IR_psth_wins(2,:,iw-1))] , wincolors(iw-1,:),...
                                    'EdgeColor','none','FaceAlpha',0.65);
                            end
                        end
                        
                        if PLOT_EX
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
                        
                        N = N+1;
                        
                        % Finish example PSTH plot
                        linkaxes(hs,'x')
                        xlim([-500 minTime])
                        
                        subplot(hs(2))
                        set(gca,'ylim',[0 add_y+1],'xtick',[])
                        
                        subplot(hs(3))
                        ylim([0 60])
                        xlabel('Time from transition (ms)')
                        ylabel('Spikes/s')
                        legend(ip,legstr)
                        
                        savedir = fullfile(fn.processed,'Transitions');
                        print_eps_kp(hfex,fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_PDCex',session,channel,clu,ib)))
%                         set(hfex,'PaperOrientation','landscape');
%                         print(hfex,'-dpdf',fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_PDCex',session,channel,clu,ib)),'-bestfit')
                        
                        end
                        
                        
                        %% Bootstrap random combinations of trials
                                                
                        IR_psth_wins = nan( numel(prevBlocks), mode(diff(t_win)), numel(t_win)-1 );
                        DistributionDiffSpikeCount = nan(numiterations,numel(t_win)-1);
                        
                        
                        for iteration = 1:numiterations
                            
                            if PLOT_EX && ~RNDsaved
                                hfexb = figure;
                                set(hfexb,'Position',figsize1,'NextPlot','add')
                                hold on
                            end
                            
                            % Permute trials to make random division in half
                            randomtrs = randperm(blocksN(ib));
                            
                            % Set up some empty vars
                            raster_x_1 = [];
                            raster_y_1 = [];
                            raster_x_2 = [];
                            raster_y_2 = [];
                            stim = nan( numel(prevBlocks), minTime-t1+1, ceil(numel(bkStart_ms)/2) );
                            psth = nan( 2, minTime-t1+1, ceil(numel(bkStart_ms)/2) );
                            
                            for ihalf = 1:2
                                
                                thesetrs = randomtrs(ihalf:2:blocksN(ib));
                                
                                % Collect FR                                 
                                for jj = 1:numel(thesetrs)
                                    psth(ihalf,:,jj) = Stream_FRsmooth( (bkStart_ms(thesetrs(jj))+t1) : (bkStart_ms(thesetrs(jj))+minTime) );
                                    stim(ihalf,:,jj) = envelope(SoundData(2, round(( (bkStart_ms(thesetrs(jj))+t1) : ((bkStart_ms(thesetrs(jj))+minTime)))  ./1000.*Info.fs_sound) ),20,'rms');
                                    sp=[]; sp = spiketimes( spiketimes>=(bkStart_ms(thesetrs(jj))+t1) & spiketimes<((bkStart_ms(thesetrs(jj))+minTime)) ) - bkStart_ms(thesetrs(jj))+1;
                                    eval( sprintf('raster_x_%i = [raster_x_%i sp];',ihalf,ihalf) )
                                    eval( sprintf('raster_y_%i = [raster_y_%i jj*ones(1,numel(sp))];',ihalf,ihalf) )
                                end
                                
                                for iw = 2:numel(t_win)
                                    % Save the psth data for each time window
                                    IR_psth_wins(ihalf,:,iw-1) = mean(psth(ihalf,(t_win(iw-1):t_win(iw)-1)-t1,:),3,'omitnan');
                                    
                                end %iw
                                
                                
                            end %ihalf
                            
                            
                            % Calculate the difference between the PSTHs for
                            % each time window
                            for iw = 2:numel(t_win)
                                DistributionDiffSpikeCount(iteration,iw-1) = sum(abs( diff(IR_psth_wins(:,:,iw-1)./1000,1) ));
                                
                                if PLOT_EX && ~RNDsaved
                                % Plot the space between PSTHs
                                hs(3)=subplot(5,1,4:5); hold on
                                fill( [t_win(iw-1)+1:t_win(iw) t_win(iw):-1:t_win(iw-1)+1] , [IR_psth_wins(1,:,iw-1) fliplr(IR_psth_wins(2,:,iw-1))] , wincolors(iw-1,:),...
                                    'EdgeColor','none','FaceAlpha',0.65);
                                end
                            end
                            
                            if PLOT_EX && ~RNDsaved
                                
                                % Now plot each PSTH
                                add_y=0;
                                for ihalf = 1:2
                                    thesetrs = randomtrs(ihalf:2:blocksN(ib));
                                    
                                    subplot(hs(3)); hold on; box off
                                    ip(ihalf)=plot(t1:minTime,mean(psth(ihalf,:,:),3,'omitnan'),'Color',ihalf.*[0.25 0.25 0.25],'LineWidth',4);
                                    hold on
                                    legstr{ihalf} = ['Random half of trials ' num2str(ihalf)];
                                    
                                    hs(1)=subplot(5,1,1); hold on; box off
                                    plot(t1:minTime,mean(stim(ihalf,:,:),3,'omitnan'),'Color',ihalf.*[0.25 0.25 0.25],'LineWidth',4)
                                    set(gca,'xtick',[],'ytick',[])
                                    
                                    hs(2)=subplot(5,1,2:3); hold on; box off
                                    plot(eval(['raster_x_' num2str(ihalf)]),eval(['raster_y_' num2str(ihalf)])+add_y,...
                                        '.','MarkerSize',15,'Color',ihalf.*[0.25 0.25 0.25])
                                    add_y = add_y + max(eval(['raster_y_' num2str(ihalf)]));
                                    
                                end
                                
                                % Finish example PSTH plot
                                linkaxes(hs,'x')
                                xlim([-500 minTime])
                                
                                subplot(hs(2))
                                set(gca,'ylim',[0 add_y+1],'xtick',[])
                                
                                subplot(hs(3))
                                ylim([0 60])
                                xlabel('Time from transition (ms)')
                                ylabel('Spikes/s')
                                legend(ip,legstr)
                                
                                
                                % if a more random combo of trials, ask to save
                                if abs( sum(ismember(thesetrs,ii)) - numel(ii)/2 )<1 && ~RNDsaved
                                    UserChoice = input('Do you want to save this plot? y/n  ::  ','s');
                                    switch UserChoice
                                        case 'y'
%                                             set(hfexb,'PaperOrientation','landscape');
%                                             print(hfexb,'-dpdf',fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_RNDex',session,channel,clu,ib)),'-bestfit')
%                                             set(hfexb,'PaperPositionMode','auto');
%                                             print(hfexb,'-painters','-depsc', fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_RNDex',session,channel,clu,ib)))
                                            savedir = fullfile(fn.processed,'Transitions');
                                            print_eps_kp(hfexb,fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_RNDex',session,channel,clu,ib)),1)
                                            
                                            RNDsaved = 1;
                                            
                                        case 'n'
                                    end
                                end
                                
                            end %if PLOT_EX
                            
                        end %iteration
                        
                        
                        %% Plot distribution of n spikes difference for bootstrapped trial combos
                        
                        if PLOT_Distr
                            
                            thisWin = 1;
                            
                            xhist=linspace(floor(min(DistributionDiffSpikeCount(:,thisWin))),ceil(max(DistributionDiffSpikeCount(:,thisWin))),100);
                            bootstrapDistr=hist(DistributionDiffSpikeCount(:,thisWin),xhist);
                            
                            hfd=figure;
                            bar(xhist,bootstrapDistr/numiterations,1,'EdgeColor',[0.7 0.7 0.7],'FaceColor',wincolors(thisWin,:))
                            hold on
                            plot(TestDiffSpikeCount(1,thisWin),0.05,'+','Color',mean(pbcolors,1),'LineWidth',2,'MarkerSize',20)
                            xlim([0 xhist(end)+1])
                            ylim([0 0.07])
                            xlabel('n spikes difference between PSTHs from bootstrapped trials')
                            ylabel('Proportion')
                            title(sprintf('IRtransition -- spike difference distribution\n%s ch%i clu%i bk%i, %i-%i ms',session,channel,clu,ib,t_win(thisWin)-1,t_win(thisWin+1)-1))
                            
                            savedir = fullfile(fn.processed,'Transitions');
                            set(hfd,'PaperOrientation','landscape');
                            print(hfd,'-dpdf',fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_distrb_%i-%ims',session,channel,clu,ib,t_win(thisWin)-1,t_win(thisWin+1)-1)),'-bestfit')
                            print(hfd,fullfile(savedir,sprintf('IRtransition_%s_%i_%i_bk%i_distrb_%i-%ims',session,channel,clu,ib,t_win(thisWin)-1,t_win(thisWin+1)-1)),'-depsc','-tiff')
                            
                        end
                        
                        
                        %% z score
                        
                        % Calculate how far the Test spike counts are from
                        % the bootstrapped distributions
                        zScores = (TestDiffSpikeCount - mean(DistributionDiffSpikeCount,1)) ./ std(DistributionDiffSpikeCount,1);
                        
                        
                        
                        %% Save result to table
                        
                        if ~any(strncmpi(Zdata.Properties.VariableNames,'zscore',6))
                            for iw = 2:numel(t_win)
                                Zdata.(sprintf('zScore_%i',t_win(iw)-1)) = 0;
                            end
                        end
                        
                        Zdata_addrow = {subject session channel clu ib [spl lpn]};
                        for iw = 2:numel(t_win)
                            Zdata_addrow{end+1} = zScores(iw-1);
                        end
                        
                        Zdata = [Zdata; Zdata_addrow];
                        
                        
                        clear zScores
                        
                        
                    end %ib
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end %session


end %subject

%%

if nargin>0
    return
end

% Save data table
Zdata(1,:) = [];
savedir = fullfile(fn.processed,'Transitions');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
writetable(Zdata,fullfile(savedir,'IRtransitionsCDF_SU'));


% Plot the result
hf=figure;
plot([2.5 2.5],[0 1],'Color',[0.5 0.5 0.5])
hold on
plot([-2.5 -2.5],[0 1],'Color',[0.5 0.5 0.5])
for iw = numel(t_win):-1:2
    [p2,x2]=ecdf(Zdata.(sprintf('zScore_%i',t_win(iw)-1)));
    ip(iw-1)=plot(x2,p2,'LineWidth',5,'Color',wincolors(iw-1,:));
    legtext{iw-1} = sprintf('%i-%i ms',t_win(iw)-1-mode(diff(t_win)),t_win(iw)-1);
end

% title(subject)
xlabel('z-score of spike count difference when trials split by prev pdc rate')
ylabel('Probability')
legend(ip,legtext,'Location','southeast','Interpreter','none')
text(6,0.5,sprintf('N = %i',size(Zdata,1)))


% Save figure

set(hf,'PaperOrientation','landscape');
print(hf,'-dpdf',fullfile(savedir,'IRtransitionsCDF_SU'),'-bestfit')
print(hf,fullfile(savedir,'IRtransitionsCDF_SU'),'-depsc','-tiff')
print(hf,fullfile(savedir,'IRtransitionsCDF_SU'),'-dsvg')

keyboard


end %function




