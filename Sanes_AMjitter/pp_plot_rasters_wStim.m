function raster = pp_plot_rasters_wStim(subject, session, channel, clu)
%
%  pp_plot_rasters(subject, session, channel, clu)  
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%
%  KP, 2016-04; last updated 2017-01


% set(0,'DefaultAxesFontSize',6)
% set(0,'defaulttextfontsize',50)

% IF SAVING PDF FILES
figFontSize      = 16;
rasterMarkerSize = 10;
rasterLineWidth  = 0.5;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;

set(0,'DefaultTextInterpreter','none')


%%%
ymaxval = 100;
%%%

% Load data files

datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
fprintf('loading data...\n')
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(datadir,subject,filename));
filename = sprintf( '%s_sess-%s_Info'  ,subject,session); load(fullfile(datadir,subject,filename));
filename = sprintf( '%s_sess-%s_Stim'  ,subject,session); load(fullfile(datadir,subject,filename));


% Manually select which flagged trials contain disruptive artifact
if ~any(strcmp(fieldnames(Info.artifact_trs),'manual'))
    addpath('helpers')
    
    Info = manually_mark_trials(Info,session,Stim);
    
    % Re-save Info structure
    filename = sprintf( '%s_sess-%s_Info',subject,session);
    save(fullfile(datadir,subject,filename),'Info','-v7.3');
    
end

%%

% GET STIM INFO

    % Get unique stimuli from this block
%     thisblk = raster(ks);
%     [ust, uIDs, stID] = unique([thisblk.fileIDs]);
    
% Get AM jitter stimuli
jStim = Stim(~isnan([Stim.fileID]'));
jtrs = find(~isnan([Stim.fileID]'));
Par_matrix = [ [jStim.block]'...
               [jStim.fileID]'...
               [jStim.dB]'...
               [jStim.behaving]'...
               [jStim.AMdepth]'...
               [jStim.HP]'...
               [jStim.LP]'];
[unique_stim, unique_IDs, StimID] = unique(Par_matrix,'rows','sorted');

% Make raster struct
raster = struct(); 
for ks = 1:size(unique_stim,1)
    
    switch unique_stim(ks,4)
        case 0, BS = 'P';
        case 1, BS = 'A';
        case 2, BS = 'D';
    end
    
    ks_trs = []; ks_trs=jtrs(StimID==StimID(unique_IDs(ks)));
    
    % Find any trials quit early
    if numel(unique([Stim(ks_trs).stimDur]))>1
        rmtr = find( [Stim(ks_trs).stimDur] ~= mode([Stim(ks_trs).stimDur]));
%         keyboard
        ks_trs(rmtr) = [];
    end
    
    raster(ks).tr_idx   = ks_trs;
    raster(ks).block    = unique_stim(ks,1);
    raster(ks).fileIDs  = unique_stim(ks,2);
    raster(ks).behaving = BS;
    raster(ks).dB       = unique_stim(ks,3);
    raster(ks).AMdepth  = unique_stim(ks,5);
    raster(ks).HP       = unique_stim(ks,6);
    raster(ks).LP       = unique_stim(ks,7);
    raster(ks).AMrate   = 4;
    raster(ks).jitter   = strtok(extractAfter(Stim(ks_trs(1)).stimfile, '4Hz_'),'.');
    raster(ks).AMonset  = round(mean([Stim(ks_trs).onsetAM] - [Stim(ks_trs).onset]));
    raster(ks).stimDur  = Stim(ks_trs(1)).stimDur;
    raster(ks).stimfn   = Stim(ks_trs(1)).stimfile;
    raster(ks).stim_str = sprintf('ch %s unit %s  |  jitter %s\n  AM %iHz  |  depth: %i pct  |  %2.3g dBSPL  |  noise band: %i-%i Hz  |  %s',...
        num2str(channel), num2str(clu), raster(ks).jitter, ...
        raster(ks).AMrate, round(raster(ks).AMdepth*100), raster(ks).dB,...
        raster(ks).HP, raster(ks).LP, raster(ks).behaving );
end

%%
% GET STIM WAVEFORM

Wave = ap_stimplotting(subject,raster);

Wcolor = [0.85 0.85 0.85];

%%

% GET SPIKE TIMES
spikes = Spikes.sorted(channel);
unit_in = find(spikes.assigns==clu);
spiketimes = round(spikes.spiketimes(unit_in) * 1000);  %ms
spiketrials = spikes.trials(unit_in); 

if isempty(spiketimes)
    error('no spike events found for this clu')
elseif spikes.labels(spikes.labels(:,1)==clu,2) == 4
    warning('  this clu is labeled as noise. are you sure you want to plot?')
    keyboard
end

% Set up figure
fprintf('plotting.\n')
nSubPlots = 2;
hS = zeros(numel(raster),nSubPlots);

for ks = 1:numel(raster)
    
%     if ~(raster(ks).block==80 && raster(ks).AMdepth==0.75 && strcmp(raster(ks).jitter,'100_rnd-1'))
%         continue
%     else
%         keyboard
%     end
    
    if numel(raster(ks).tr_idx)<3
        continue
    end
       
    % Set current figure/subplot handles
    hF(ks) = figure; hold on
    scrsz = get(0,'ScreenSize');
    set(hF(ks),'Position',[1 scrsz(4) scrsz(3) scrsz(4)],...
        'Nextplot','add');
    for isp = 1:nSubPlots
        hS(ks,isp)=subplot(nSubPlots,1,isp);
        set(hS(ks,isp),'Nextplot','add');
    end
    
    
    % Set up raster/histo plot parameters
    bin    = 20;    %ms
    t_beg  = -399;  %ms
    try
        t_end  = floor( min(400+raster(ks).stimDur, Info.t_win_ms(2)) /bin)*bin;  %ms  (prefer this: max([raster.stimDur])+400  but tricky to reshape)
    catch
        keyboard
    end
    nt     = t_end - t_beg +1;  %each entry 1 ms
    
    smooth.wsize = round(nt/200);   %window size for gaussian smoothing of histo for plotting
    smooth.cutoff = 20;   %cutoff for the gaussian smoothing
    smooth.stdev = Info.fs/(2*pi*smooth.cutoff); %std for gaussian smoothing of histo for plotting


    raster(ks).window_ms = [t_beg t_end];
    
    
    % Get trials for this stim and remove flagged trials
    tr_this_stim = raster(ks).tr_idx;
    
    FLAGGED = [];
    FLAGGED = Info.artifact_trs(channel).manual;
    
    
    [~,ibt,~] = intersect(tr_this_stim,FLAGGED);
    tr_this_stim(ibt) = [];
    
    
    % Get spiketimes and create raster
    raster_x=[];  raster_y=[];  hist_raw=zeros(1,nt);
    for it = 1:numel(tr_this_stim)
        
        sp=[];  spk_in=[];
        spk_in = find(spiketrials==tr_this_stim(it));
        sp = spiketimes(spk_in) + ones(size(spiketimes(spk_in)))*(Info.t_win_ms(1)-1); %ms, rel to t0
        sp = sp( sp>=t_beg & sp<= t_end );
        
        hist_raw(it, sp-t_beg+1) = 1;
        raster_x = [raster_x sp];
        raster_y = [raster_y it .* ones(1,numel(sp))];
    end
    
    raster(ks).x = raster_x;
    raster(ks).y = raster_y;
    
    hist_raw = sum(hist_raw,1) / it;
    hist_bin = sum(reshape(hist_raw, bin, nt/bin),1)/(bin/1000);
    hist_smooth = smoothts(hist_bin,'g', smooth.wsize, smooth.stdev);
    
    
    
    % Plot this stimulus
        
    % raster
    subplot(hS(ks,1)); hold on
%     fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*(1+it) + (1+it)/2 ,...
%         Wcolor,'EdgeColor','none')
    plot([0 0],[-30 30],'k:', 'LineWidth', 2)
    plot([raster(ks).AMonset raster(ks).AMonset],[-30 30],'k:', 'LineWidth', 2)
    plot([raster(ks).stimDur raster(ks).stimDur],[-30 30],'k:', 'LineWidth', 2)
    % long ticks
    plot(  raster_x  ,  raster_y  , 'k+','MarkerSize',rasterMarkerSize, 'LineWidth', rasterLineWidth)
    for ii=1:it
        plot([t_beg t_end], [ii ii], 'k', 'LineWidth', rasterLineWidth)
    end
    axis tight
    set(gca, 'XLim', [t_beg t_end], 'XTick',[], 'YLim', ([0 1+it]))
    ylabel('Trials')
    set(gca,'FontSize',figFontSize)
    hold off
        
    
    % psth
    subplot(hS(ks,2)); hold on
    fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*ymaxval + (ymaxval)/2 ,...
        Wcolor,'EdgeColor','none')
    plot( t_beg:bin:t_end  , hist_bin , 'k', 'LineWidth', 2)
    plot([0 0],[0 ymaxval],'k:', 'LineWidth', 2)
    plot([raster(ks).AMonset raster(ks).AMonset],[0 ymaxval],'k:', 'LineWidth', 2)
    plot([raster(ks).stimDur raster(ks).stimDur],[0 ymaxval],'k:', 'LineWidth', 2)
    set(gca, 'XLim', [t_beg t_end])
    xlabel( 'Time (ms)')
    ylabel('Spikes/sec')
    ylim([0 ymaxval])
    set(gca,'FontSize',figFontSize)
    hold off
    
    % Set font sizes
    set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)
    
    % Add title
    t=suptitle(raster(ks).stim_str);
    t.FontSize=0.75*figFontSize;
    
    
    % SAVE FIGURE
    
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
    if raster_y<6
        savedirR = fullfile(savedir,subject,'^rasters',session,['ch' num2str(channel)],'tested');
    else
        savedirR = fullfile(savedir,subject,'^rasters',session,['ch' num2str(channel)]);
    end
    if ~exist(savedirR,'dir')
        mkdir(savedirR)
    end
    savename = sprintf('%s_%s_ch%i_clu%i_AM%iHz_jitter%s_%idpth_%idB_%i-%i_%s_blk%i',subject,session,channel,clu,...
        raster(ks).AMrate, raster(ks).jitter, round(raster(ks).AMdepth*100),...
        raster(ks).dB, raster(ks).HP, raster(ks).LP, raster(ks).behaving, raster(ks).block);
%     print(hF(ks),fullfile(savedir,subject,'^rasters',session,savename),'-depsc','-tiff')
    set(gcf,'PaperOrientation','landscape');
    print(hF(ks),'-dpdf',fullfile(savedirR,savename),'-bestfit')
    
    if ks==60
        close all
    end
    
    
    % Calculate some basic response properties
    raster(ks).nSpk = sum(raster_x > raster(ks).AMonset & raster_x < raster(ks).stimDur)...
        / it / ((raster(ks).stimDur-raster(ks).AMonset)/1000); 
        %count spikes during AM, divide by n trials, divide by seconds duration of AM stim
    
    
end

% SAVE RASTER STRUCT

% savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
% save(fullfile(savedir,subject,savename),'raster','-v7.3')









end

