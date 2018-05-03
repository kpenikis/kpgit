function ap_mph_v2(subject, session, channels, clus)
%
%  pp_mph(subject, session, channel, clu)
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%


% done with plotting; now change data collection window to fixed time
%   instead of one period
% response metrics: VS, RS, mean phase, nspks, TS?


%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10; 
%!!!!!!!!!!!!!!!!!

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];
figsize2 = [1 scrsz(4) scrsz(3)/2 scrsz(4)];


histbinsize = 20;
anbinsize   = 10;
smthbinsize = 20;


colors = hsv(6);
colors = [colors; 0.5.*hsv(4)];


%%  SESSIONS

% Get list of sessions to check for sorted data

fn = set_paths_directories(subject);
SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));

try
    
Sessions = [];
for ifn = 1:numel(SpkFns)
    if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
        Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
    end
end

catch
    keyboard
end

% Step through each session
for sess = Sessions'
    
session = char(sess);


%%
% Load data files
fn = set_paths_directories(subject,session);
fprintf('Session %s:\n',session)
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


% STEP THROUGH EACH CHANNEL
for channel = channels
    
    % Artifact for this channel
    ArtifactFlag = Info.artifact(channel).SDsamples;
    
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if all(spikes.labels(:,2)==1)
            continue
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
        end
    end
    
    
    % STEP THROUGH EACH CLU
    
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
        
        fprintf('ch %i clu %i - ',channel,clu)
        if numel(spiketimes) < round(FRcutoff*length(SoundData)/Info.fs_sound)
            fprintf('skip\n')
            continue
        else
            fprintf('plotting...\n')
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
        % max of first spiketime or 5 seconds before unmodulated sound came on
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
                    
                    % Get this unit/stimulus combo's N clean blocks (trials)
                    blocks_N = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    if min(blocks_N)<minTrs
                        continue
                    end
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Create empty vectors for data
                    
                    MPH = struct();
                    for ir = 1:numel(AMrates)
                        MPH(ir).pdc.raster  = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                        MPH(ir).IR7.raster  = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                        MPH(ir).IR8.raster  = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                        MPH(ir).IR9.raster  = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                        MPH(ir).IR10.raster = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                    end
                    
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % First collect data from IR periods
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for ib = 7:10
                        
                        % Get onsets for this block type
                        [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                            ib,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                        
                        %
                        trials = randperm(numel(bkStart_ms));
                        
                        for it = 1:numel(trials)
                            
                            if it>min(blocks_N), continue, end
                            
                            newRates = -1+bkStart_samps(trials(it)) + find(diff( [SoundData(1, (bkStart_samps(trials(it))-1):(bkStop_samps(trials(it))) ) 0] ));
                            if numel(newRates)~=13
                                keyboard
                            end
                            
                            for ir = 1:(numel(newRates)-1)
                                this_rate = SoundData(1, newRates(ir));
                                t = round( [newRates(ir) newRates(ir+1)-1] / Info.fs_sound*1000 );
                                sp=[]; sp = spiketimes( spiketimes>=t(1) & spiketimes<=t(2) ) - t(1) +1;
                                
                                % Determine whether this period was in the
                                % first or second sequence of the IR block
                                which_seq = 1+ceil( (newRates(ir)-bkStart_samps(trials(it))) /Info.fs_sound - sum(1./AMrates)+0.005 );
                                if (which_seq<1 || which_seq>2)
                                    keyboard
                                end
                                
                                MPH(AMrates==this_rate).(sprintf('IR%i',ib)).pdtime(min(blocks_N)*(which_seq-1)+it,1)  = which_seq;
                                MPH(AMrates==this_rate).(sprintf('IR%i',ib)).pdtime(min(blocks_N)*(which_seq-1)+it,2)  = 1000*(newRates(ir)-bkStart_samps(trials(it))) /Info.fs_sound;
                                MPH(AMrates==this_rate).(sprintf('IR%i',ib)).pdtime(min(blocks_N)*(which_seq-1)+it,3)  = SoundData(8, bkStart_samps(trials(it))-1 );
                                
                                % Put spikes into corresponding raster
                                MPH(AMrates==this_rate).(sprintf('IR%i',ib)).raster(min(blocks_N)*(which_seq-1)+it,sp) = 1;
                                
                                
                            end
                            
                        end %it
                        
                    end %ib
                    
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Now collect data from periodic blocks
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for ib = 1:6
                        
                        % Get onsets for this block type
                        [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                            ib,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                        
                        trials = randperm(numel(bkStart_ms));
                        
                        for it = 1:numel(trials)
                            if it>min(blocks_N), continue, end
                            
                            % Estimate samples of period restarts
                            % (phase=0)
                            pd_starts = 0 : round(Info.fs_sound/AMrates(ib)) : (bkStop_samps(trials(it))-bkStart_samps(trials(it)));
                            
                            % Get closest period to start before/after these times (in seconds)
                            which_pd = [0 1];%[0 0.5 1 1.5];
                            
                            for ipd = 1:length(which_pd)
                                [~,pdc_pd] = min(abs(pd_starts-which_pd(ipd)*Info.fs_sound));
                                pd_start_ms = round( ( pd_starts(pdc_pd) + bkStart_samps(trials(it)) ) /Info.fs_sound *1000 );
                                %                                 eval(sprintf('MPH_pdc_%i(it, spiketimes( spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) - pd_start_ms +1 ) = 1;;',AMrates(ib) ))
                                
                                MPH(ib).pdc.pdtime(min(blocks_N)*(ipd-1)+it,1)  = ipd;
                                MPH(ib).pdc.pdtime(min(blocks_N)*(ipd-1)+it,2)  = round( pd_starts(pdc_pd)/Info.fs_sound *1000 );
                                MPH(ib).pdc.pdtime(min(blocks_N)*(ipd-1)+it,3)  = SoundData(8, bkStart_samps(trials(it))-1 );
                                
                                % Get spikes and put into corresponding raster
                                sp=[]; sp = spiketimes(  spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) - pd_start_ms +1 ;
                                MPH(ib).pdc.raster(min(blocks_N)*(ipd-1)+it,sp) = 1;
                                
                            end
                            
                        end %it
                    end %ib
                    
                    
                    
                    
                    %% Plots
                    
                    for ir = 2%1:numel(AMrates)
                        
                        % Make figure
                        hf(ir) = figure;
                        set(hf(ir),'Position',figsize2,'NextPlot','add')
                        
                        
                        % Make string for title
                        ratestr = [num2str(AMrates(ir)) 'Hz'];
                        unType = {'' 'SU' 'MU'};
                        titlestr1 = sprintf('MPHs for %s\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz',...
                            ratestr,subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,SoundData(5,bkStart_samps(1)),lpn);
                        
                        % Get data for this rate from each context
                        % (pdc, IR 7:10, first/middle)
                        contextcolors = [colors(ir,:); colors(7:10,:)];
                        nrows = 0;
                        max_count = 0;
                        
                        context = fieldnames(MPH);
                        for icontext = context'
                            for ii = 1:numel(unique(MPH(ir).pdc.pdtime(:,1)))
                                
                                these_rows = find(MPH(ir).(icontext{:}).pdtime(:,1)==ii)';
                                
                                % Get raster data
                                raster_x = []; raster_y = [];
                                for it = 1:numel(these_rows)
                                    raster_x = [raster_x find( MPH(ir).(icontext{:}).raster(these_rows(it),:) )];
                                    raster_y = [raster_y repmat(nrows+it,1,sum( MPH(ir).(icontext{:}).raster(these_rows(it),:) ))];
                                end
                                
                                % Plot rasters (left early, right middle)
                                ax_r(ii)=subplot(6,2,ii);
                                hold on
                                hL=plot(raster_x, raster_y,'+','Color',contextcolors(strcmp(icontext,context),:),...
                                    'MarkerSize',3,'LineWidth',0.5);
                                % If first column, set axis labels
                                if ii==1
                                    xlabel('time (ms)')
                                    ylabel([num2str(min(blocks_N)) ' trs each'])
                                end
%                                 % Set marker transparency
%                                 if ~isempty(hL)
%                                     hMarkers=[];
%                                     while isempty(hMarkers)
%                                         hMarkers = hL.MarkerHandle;
%                                         pause(0.1)
%                                     end
%                                     hMarkers.EdgeColorData(4) = uint8(255 * 0.4 );
%                                 end

                                
                                
                                % If IR context, separate by preceding block
                                bardata = [];
                                if strcmp(icontext{:}(1:2),'IR') && numel(unique(MPH(ir).(icontext{:}).pdtime(:,3)))==2
                                    bardata = [];
                                    for ip = unique(MPH(ir).(icontext{:}).pdtime(:,3))'
                                        bardata = [bardata; hist(raster_x(ismember(raster_y-nrows,find(MPH(ir).(icontext{:}).pdtime(:,3)==ip))),linspace(0,1000/AMrates(ir),52))];
%                                         raster_x(ismember(raster_y,find(MPH(ir).(icontext{:}).pdtime(:,3)==ip)))
                                    end
                                else
                                    bardata = hist(raster_x,linspace(0,1000/AMrates(ir),52));
                                end
                                
                                % Plot folded MPH
                                
                                sp_idx = 2*find(strcmp(icontext,context)) + ii;
                                ax_m(sp_idx-2) = subplot(6,2,sp_idx);
                                hold on
                                hBar = bar(bardata',1,'stacked');
                                barcols = [contextcolors(strcmp(icontext,context),:); 0 0 0];
                                for ip = 1:size(bardata,1)
                                    set(hBar(ip),'FaceColor',barcols(ip,:),'EdgeColor',barcols(ip,:));
                                end
%                                     'FaceColor',contextcolors(strcmp(icontext,context),:),...
%                                     'EdgeColor',contextcolors(strcmp(icontext,context),:))
                                box off
                                title([icontext{:} ', ' num2str(round(mode(MPH(ir).(icontext{:}).pdtime(these_rows,2)))) 'ms'])
                                
                                % Keep track of largest bin spk count to
                                % set ymax
                                if max(hist(raster_x,linspace(0,1000/AMrates(ir),52))) > max_count
                                    max_count = max(hist(raster_x,linspace(0,1000/AMrates(ir),52)));
                                end
                                
                            end %ii
                            
                            % Add to count of total number of rows
                            nrows = nrows+min(blocks_N);
                            
                        end %context
                        
                        % Link subplot axes and set limits
                        linkaxes(ax_r,'xy')
                        linkaxes(ax_m,'xy')
                        set(ax_r,'xlim',[0 1000/AMrates(ir)],'ylim',[0 nrows])
                        set(ax_m,'xlim',[0.5 52.5],'ylim',[0 max_count+1])
                        
                        % Adjust axis labels
                        set(ax_m,'xticklabel',[])
                        set(ax_m(2:end),'yticklabel',[])
                        set(ax_r,'yticklabel',[])
                        
                        % Add overall title
                        suptitle(titlestr1)
                        
                        
                        % SAVE FIGURE
                        
%                         savedir = fullfile(fn.processed,subject,'^an_plots','MPH');
%                         if ~exist(savedir,'dir')
%                             mkdir(savedir)
%                         end
%                         
%                         savename = sprintf('%s_%s_ch%i_clu%i_%s_%idB_LP%ihz_MPH_%iHz',...
%                             subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn,AMrates(ir));
%                         set(hf(ir),'PaperOrientation','landscape');
%                         print(hf(ir),'-dpdf',fullfile(savedir,savename),'-bestfit')
%                         
                        
                    end %ir
                    
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end %session

end %function




