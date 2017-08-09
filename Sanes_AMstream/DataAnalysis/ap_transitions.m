function ap_transitions(subject,session,channels,clus)
%
%  ap_zscore_plots(subject, session, channel, clu)
%    Similar to ap_plot_rasters, but converts FR to a zscore.
%
%  KP, 2017-08
 

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];

histbinsize = 20;
anbinsize   = 50;
smthbinsize = 50;

colors = hsv(6);
colors = [colors; 0.5.*hsv(4)];

unType = {'' 'SU' 'MU'};

N = 0;


%%  SESSIONS

% Get list of sessions to check for sorted data

% fn = set_paths_directories(subject);
% SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
% 
% try
%     
% Sessions = [];
% for ifn = 1:numel(SpkFns)
%     if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
%         Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
%     end
% end
% 
% catch
%     keyboard
% end
% 
% % Step through each session
% for sess = Sessions'
%     
% session = char(sess);

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

AMrates = [2 4 8 16 32 64];

TotalDur_s = size(SoundData,2)/Info.fs_sound;


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
%             disp(' SESSION MAY NOT BE MANUALLY SORTED YET')
            continue
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
%             disp(' !! no valid clus for this channel')
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = clus'
        
        % !! Only SU for now !!
        if spikes.labels(spikes.labels(:,1)==clu,2) ~= 2
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
        elseif numel(spiketimes) < round(3*length(SoundData)/Info.fs_sound)
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
                    
%                     for ir = AMrates
%                         % Find all blocks that start with this rate, no
%                         % matter what follows
%                             % have to store each if want to add to it later,
%                             % from other blocks (if analyzing just a given
%                             % period, as opposed to a specified duration
%                             % from the onset of a block)
%                     end

                    % Get this unit/stimulus combo's N clean blocks (trials)
                    blocks_N = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    t1 = -999;
                    t2 = 1000;

                    % Go through each stiulus
                    for ib = 1:numel(Info.blockKey)
                        
                        %if unmodulated or silent, skip for now
                        if ib>=11, continue, end
                        
                        % Make figure
                        hf = figure;
                        set(hf,'Position',figsize1,'NextPlot','add')
                        hold on
                        
                        % Get this block start times
                        [bkStart_samps,bkStop_samps,...
                            bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                                                ib,spl,lpn,amd,ArtifactFlag,fs);
                        
                        %view preceeding blocks:
%                         hpb = hist(SoundData(8,bkStart_samps-1),0:12);
%                         hpr = hist(SoundData(1,bkStart_samps-1),[0 3 5 9 17 33 65]);
                        
                        % Get list of blocks that came before this one
                        prevBlocks = sort(unique(SoundData(8,bkStart_samps-1)));
                        
                        % Set up some empty vars
                        legstr = cell(1,numel(prevBlocks)); clear ip
                        stim = nan(numel(bkStart_samps),max(bkStop_samps-bkStart_samps));
                        
                        for pb = prevBlocks
                            
                            ii = find(SoundData(8,bkStart_samps-1) == pb);
                            
                            % Collect FR for specified transitions
                            psth_smooth = nan(numel(bkStart_ms(ii)),max(bkStop_ms(ii)-bkStart_ms(ii)));
                            
                            for jj = 1:numel(ii)
                                psth_smooth(jj,1:(bkStop_ms(ii(jj))-bkStart_ms(ii(jj))-t1))    = Stream_FRsmooth((bkStart_ms(ii(jj))+t1):(bkStop_ms(ii(jj))-1));
                                stim(ii(jj),1:floor(bkStop_samps(ii(jj))-bkStart_samps(ii(jj))-t1/1000*fs)) = envelope(SoundData(2, ceil(bkStart_samps(ii(jj))+t1/1000*fs):(bkStop_samps(ii(jj))-1) ),50,'rms');
                            end
                            
                            % Plot data
                            subplot(4,1,2:4); hold on
                            fill( [t1:(t2-1) (t2-1):-1:t1], [(mean(psth_smooth(:,1:(t2-t1)),1) - std(psth_smooth(:,1:(t2-t1)),1)/sqrt(jj)) fliplr(mean(psth_smooth(:,1:(t2-t1)),1) + std(psth_smooth(:,1:(t2-t1)),1)/sqrt(jj))],...
                                colors(pb,:),'EdgeColor','none','FaceAlpha',0.6)
                            ip(pb==prevBlocks) = plot(t1:(t2-1),mean(psth_smooth(:,1:(t2-t1)),1),'Color',colors(pb,:),'LineWidth',jj/2);
                            legstr{pb==prevBlocks} = [Info.blockKey{pb} ', n=' num2str(jj)];
                            
                        end
                        
                        
                        % Plot formatting
                        subplot(4,1,2:4); hold on
                        if ib<7
                            blockstr = [Info.blockKey{ib} 'Hz'];
                        else
                            blockstr = Info.blockKey{ib};
                        end
                        legend(ip,legstr)
                        xlim([t1 t2])
                        xlabel('time from transition (ms)')
                        ylabel('spikes/sec')
                        
                        % Plot stimulus amplitude above
                        subplot(4,1,1)
                        plot(stim(:,round([1:(t2-t1)]/1000*fs))','Color',colors(ib,:),'LineWidth',2);
                        xlim([1 (t2-t1)])
                        set(gca,'xtick',[],'ytick',[])
                        
                        % and a title
                        titlestr1 = sprintf('Block %s\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz',...
                            blockstr,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn);
                        title(titlestr1)
                        
                        
                        % Save the figure
                        
                        
                        
                    end %ib
                    
                    
                    %% 
                    % Perform comparison of each sequence, 
                    % when it comes early or late in the block
                                        
                    for x = {'IR_A' 'IR_B' 'IR_C' 'IR_D'}
                        
                        clear ip
                        legstr = cell(1,2);
                        
                        hf=figure;
                        set(hf,'Position',figsize1,'NextPlot','add')
                        hold on
                        
                        % Which blocks contain this sequence?
                        blocks = find(~cellfun(@isempty,strfind(Info.blockKey,x{:})));
                                                
                        for ib = blocks
                            
                            % Get this block start times
                            [bkStart_samps,bkStop_samps,...
                                bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                                                ib,spl,lpn,amd,ArtifactFlag,fs);
                            
                            % Did the sequence come first or second in this block?
                            switch strfind(Info.blockKey{ib},x{:})
                                case 1
                                    pos = 'first';
                                    linecol = [0.8 0.6 0.2];
                                case 6
                                    pos = 'second';
                                    linecol = [0.8 0.4 0.4];
                            end
                            
                            stim = nan(numel(bkStart_ms),ceil(max(bkStop_samps-bkStart_samps)/2));
                            psth_smooth = nan(numel(bkStart_ms),1000);
                            for itr = 1:numel(bkStart_ms)
                                                                
                                % Find time that second sequence begins
                                pdOnsets = find(diff( SoundData(1,(bkStart_samps(itr)-1):bkStop_samps(itr)) ));
                                samp_2ndSeqOnset = pdOnsets(7)+bkStart_samps(itr);
                                
                                switch pos
                                    case 'first'
                                        seqsamps = bkStart_samps(itr):(samp_2ndSeqOnset-1);
                                    case 'second'
                                        seqsamps = samp_2ndSeqOnset:bkStop_samps(itr);
                                end
                                
                                seqms = round([seqsamps(1) seqsamps(end)]/fs*1000);
                                
                                psth_smooth(itr,1:diff(seqms)) = Stream_FRsmooth(seqms(1):(seqms(2)-1));
                                stim(itr,1:(seqsamps(end)-seqsamps(1)+1)) = envelope(SoundData(2,seqsamps),50,'rms');
                                
                            end %itr
                            
                            % Plot this one
                            subplot(4,1,2:4); hold on
                            fill( [1:diff(seqms) diff(seqms):-1:1], [(mean(psth_smooth(:,1:diff(seqms)),1,'omitnan') - std(psth_smooth(:,1:diff(seqms)),1,'omitnan')/sqrt(itr)) fliplr(mean(psth_smooth(:,1:diff(seqms)),1,'omitnan') + std(psth_smooth(:,1:diff(seqms)),1,'omitnan')/sqrt(itr))],...
                                linecol,'EdgeColor','none','FaceAlpha',0.6)
                            ip(ib==blocks)=plot(1:diff(seqms),mean(psth_smooth(:,1:diff(seqms)),1,'omitnan'),'Color',linecol,'LineWidth',2);
                            legstr{ib==blocks} = [pos ', n=' num2str(itr)];
                            
                            % Plot stimulus amplitude above
                            subplot(4,1,1); hold on
                            plot(stim','Color',linecol,'LineWidth',0.75);
                            
                        end %ib
                        
                        % Add title
                        title(x{:})
                        xlim([1 size(stim,2)])
                        set(gca,'xtick',[],'ytick',[])
                        
                        % Finish data subplot
                        subplot(4,1,2:4); hold on
                        legend(ip,legstr)
                        xlim([1 diff(seqms)])
                        
                        
                        % Save figure
                        
                        
                        
                        
                    end %x
                    
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

return
% end % sessions


%% FINISH AND SAVE FIGURES

% Add title
switch subject
    case 'WWWf_244303'
        region = 'caudal A1';
    case 'WWWr_244300'
        region = 'DP/VP';
end
figure(hf1); title(sprintf('%s\nIR responses  |  N = %i',region,N))
figure(hf2); title(sprintf('%s\nIR responses  |  N = %i',region,N))


% Set save directory

savedir = fullfile(fn.processed,subject,'^an_plots','Population');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% Save figs

% mean of IR blocks
savename = sprintf('%s_allSU_%idB_LP%ihz_ObsVsPred-meanIR',...
    subject,spl,lpn);
set(hf1,'PaperOrientation','landscape');
print(hf1,'-dpdf',fullfile(savedir,savename),'-bestfit')

% each IR sequence plotted
savename = sprintf('%s_allSU_%idB_LP%ihz_ObsVsPred-eachIR',...
    subject,spl,lpn);
set(hf2,'PaperOrientation','landscape');
print(hf2,'-dpdf',fullfile(savedir,savename),'-bestfit')


end %function




