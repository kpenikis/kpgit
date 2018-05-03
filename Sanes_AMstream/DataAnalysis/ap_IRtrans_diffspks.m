function ap_IRtrans_diffspks(select_subject,select_session,channels,clus)
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


set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];

histbinsize = 5;
smthbinsize = 20;

colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors = [colors; 0.7.*bone(4)];


unType = {'' 'SU' 'MU'};

N = 0;



    
Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.IRblock  = 0;
Zdata.stimpars = [0 0];
Zdata.zScore_1 = 0;
Zdata.zScore_2 = 0;



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
        case 'WWWr_244300'
            region = 'DP/VP';
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
        
        %or with sliding 50 ms boxcar
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
                    
                    
                    %%
                    % Set up figure
                    fprintf('plotting ch %i clu %i...\n',channel,clu)

                    % Get this unit/stimulus combo's N clean blocks (trials)
                    [blocksN,minTime] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd,7:10);
                    
                    
                    %~~~~~~~~~~~~~~
                    t1 = -999;
                    t2 = 1;
                    t3 = 1000;
                    t4 = minTime;
                    %~~~~~~~~~~~~~~
                    
                    % Go through each stiulus
                    for ib = 7:10
                        
                        
                        % Make figure
                        hf = figure;
                        set(hf,'Position',figsize1,'NextPlot','add')
                        hold on
                        
                        if blocksN(ib)<minTrs, continue, end
                        
                        % Get this block start times
                        [bkStart_samps,bkStop_samps,...
                            bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                                                ib,spl,lpn,amd,ArtifactFlag,fs);
                        
                        % Get list of blocks that came before this one
                        prevBlocks = sort(unique(SoundData(8,bkStart_samps-1)));
                        
                        % Set up some empty vars
                        legstr = cell(1,numel(prevBlocks)); clear ip
                        stim = nan(numel(bkStart_samps),max(bkStop_samps-bkStart_samps));
                        
                        IR_psth_1 = nan(2,t3-t2+1);
                        IR_psth_2 = nan(2,t4-t3+1);
                        
                        for pb = prevBlocks
                            
                            ii = find(SoundData(8,bkStart_samps-1) == pb);
                            
                            % Collect FR for specified transitions
                            psth_smooth = nan(numel(bkStart_ms(ii)), t4-t1+1 );
                            
                            for jj = 1:numel(ii)
                                psth_smooth(jj, 1:((bkStart_ms(ii(jj))+minTime)-bkStart_ms(ii(jj))-t1) )...
                                    = Stream_FRsmooth( (bkStart_ms(ii(jj))+t1) : ((bkStart_ms(ii(jj))+minTime)-1) );
%                                 stim(ii(jj),1:floor(bkStop_samps(ii(jj))-bkStart_samps(ii(jj))-t1/1000*fs)) = envelope(SoundData(2, ceil(bkStart_samps(ii(jj))+t1/1000*fs):(bkStop_samps(ii(jj))-1) ),50,'rms');
                            end
                            
                            % Save the pstch data to compare to the other
                            IR_psth_1(pb==prevBlocks,:) = mean(psth_smooth(:,(t2-t1):(t3-t1)),1);
                            IR_psth_2(pb==prevBlocks,:) = mean(psth_smooth(:,(t3-t1):(t4-t1)),1);
                            
                            % Plot data
                            %  **Now use plot to show different between PSTHs
%                             subplot(4,1,2:4); hold on
                            ip(pb==prevBlocks) = plot(t1:(t4-1),mean(psth_smooth(:,1:(t4-t1)),1),'Color',colors(pb,:),'LineWidth',jj/4);
                            legstr{pb==prevBlocks} = [Info.blockKey{pb} ', n=' num2str(jj)];
                            
                            
                        end %pb
                        
                        
                        TestDiffSpikeCount_1 = sum(abs(diff(IR_psth_1./1000,1)));
                        TestDiffSpikeCount_2 = sum(abs(diff(IR_psth_2./1000,1)));
                        %correct for difference in total time in each epoc
                        TestDiffSpikeCount_2 = TestDiffSpikeCount_2/(t4-t3+1)*(t3-t2+1);
                        
                        
                        
                        %% Bootstrap random combinations of trials
                        
                        IR_psth_1 = nan(2,t3-t2+1);
                        IR_psth_2 = nan(2,t4-t3+1);
                        
                        %~~~~~~~~~~~~~~~~~~~~~~
                        numiterations = 1000;
                        %~~~~~~~~~~~~~~~~~~~~~~
                        
                        DistributionDiffSpikeCount_1 = nan(numiterations,1);
                        DistributionDiffSpikeCount_2 = nan(numiterations,1);
                        
                        for iteration = 1:numiterations
                            
                            randomtrs = randperm(blocksN(ib));
                            
                            for ihalf = 1:2
                                
                                try
                                thesetrs = randomtrs(ihalf:2:blocksN(ib));
                                catch
                                    keyboard
                                end
                                
                                % Collect FR for specified transitions
                                psth_smooth = nan(numel(bkStart_ms(thesetrs)),max(bkStop_ms(thesetrs)-bkStart_ms(thesetrs)));
                                
                                for jj = 1:numel(thesetrs)
                                    psth_smooth(jj,1:(bkStop_ms(thesetrs(jj))-bkStart_ms(thesetrs(jj))-t1))    = Stream_FRsmooth((bkStart_ms(thesetrs(jj))+t1):(bkStop_ms(thesetrs(jj))-1));
                                    stim(thesetrs(jj),1:floor(bkStop_samps(thesetrs(jj))-bkStart_samps(thesetrs(jj))-t1/1000*fs)) = envelope(SoundData(2, ceil(bkStart_samps(thesetrs(jj))+t1/1000*fs):(bkStop_samps(thesetrs(jj))-1) ),50,'rms');
                                end
                                
                                % Save the pstch data to compare to the other
                                IR_psth_1(ihalf,:) = mean(psth_smooth(:,(t2-t1):(t3-t1)),1);
                                IR_psth_2(ihalf,:) = mean(psth_smooth(:,(t3-t1):(t4-t1)),1);
                                
%                             subplot(4,1,2:4); hold on
%                             ip(ihalf) = plot(t1:(t4-1),mean(psth_smooth(:,1:(t4-t1)),1),'Color',colors(ihalf,:),'LineWidth',jj/4);
%                             legstr{pb==prevBlocks} = [Info.blockKey{pb} ', n=' num2str(jj)];
                            
                            end %ihalf
                            
                            
                            DistributionDiffSpikeCount_1(iteration) = sum(abs(diff(IR_psth_1./1000,1)));
                            DistributionDiffSpikeCount_2(iteration) = sum(abs(diff(IR_psth_2./1000,1)));
                            DistributionDiffSpikeCount_2(iteration) = DistributionDiffSpikeCount_2(iteration)/(t4-t3+1)*(t3-t2+1);
                            
                            
                        end %iteration
                        
                        
                        
                        %% z score
                        
                        % Calculate how far the Test spike counts are from
                        % the bootstrapped distributions
                        zScore_1 = (TestDiffSpikeCount_1 - mean(DistributionDiffSpikeCount_1)) / std(DistributionDiffSpikeCount_1);
                        zScore_2 = (TestDiffSpikeCount_2 - mean(DistributionDiffSpikeCount_2)) / std(DistributionDiffSpikeCount_2);
                        
                        if zScore_1>4
                            keyboard
                        end
                        
                        
                        %% Save result to table
                        
                        Zdata_addrow = {subject session channel clu ib [spl lpn] zScore_1 zScore_2};
                        
                        Zdata = [Zdata; Zdata_addrow];
                        
                        
                        clear zScore_1 zScore_2
                        
                        
                        %%
%                         % Plot formatting
%                         subplot(4,1,2:4); hold on
%                         if ib<7
%                             blockstr = [Info.blockKey{ib} 'Hz'];
%                         else
%                             blockstr = Info.blockKey{ib};
%                         end
%                         legend(ip,legstr)
%                         xlim([t1 t3])
%                         xlabel('time from transition (ms)')
%                         ylabel('spikes/sec')
%                         
%                         % Plot stimulus amplitude above
%                         subplot(4,1,1)
%                         plot(stim(:,round([1:(t3-t1)]/1000*fs))','Color',colors(ib,:),'LineWidth',2);
%                         xlim([1 (t3-t1)])
%                         set(gca,'xtick',[],'ytick',[])
%                         
%                         % and a title
%                         titlestr1 = sprintf('Block %s\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz',...
%                             blockstr,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn);
%                         title(titlestr1)
                        
                        
%                         %% Save figure
%                         if ib>6
%                             
%                             savedir = fullfile(fn.processed,subject,'^an_plots','Transitions');
%                             if ~exist(savedir,'dir')
%                                 mkdir(savedir)
%                             end
%                             
%                             savename = sprintf('%s_%s_ch%i_clu%i_%s_%idB_LP%ihz_%s',...
%                                 subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn,Info.blockKey{ib});
%                             set(hf,'PaperOrientation','landscape');
%                             print(hf,'-dpdf',fullfile(savedir,savename),'-bestfit')
%                             
%                         end
                        
                    end %ib
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end %session


end %subject

%%

% Plot the result

hf=figure;
plot([2.5 2.5],[0 1],'Color',[0.5 0.5 0.5])
hold on
[p2,x2]=ecdf(Zdata.zScore_2);
ip(2)=plot(x2,p2,'Color',[139 71 38]./255,'LineWidth',6);
[p1,x1]=ecdf(Zdata.zScore_1);
ip(1)=plot(x1,p1,'Color',[255 120 3]./255,'LineWidth',6);

% title(subject)
xlabel('z-score of spike count difference when trials split by prev pdc rate')
ylabel('Probability')
legend(ip,{'first 1000 ms' 'second 1000 ms'},'Location','southeast')
text(6,0.25,sprintf('N = %i (%i units)',size(Zdata,1),size(Zdata,1)/4))


% Save figure

savedir = fullfile(fn.processed,'Transitions');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
set(hf,'PaperOrientation','landscape');
print(hf,'-dpdf',fullfile(savedir,'IRtransitionsCDF_SU'),'-bestfit')


% Save data table
Zdata(1,:) = [];
writetable(Zdata,fullfile(savedir,'IRtransitionsCDF_SU'));




end %function




