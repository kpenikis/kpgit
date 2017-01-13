function ap_plot_TC_AMrate(subject, session, channel, clu)
%
%  pp_plot_rasters(subject, session, channel, clu)  
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%
%  KP, 2016-04; last updated 2016-12


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

datadir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
fprintf('loading data...\n')
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(datadir,subject,filename));
filename = sprintf( '%s_sess-%s_Info'  ,subject,session); load(fullfile(datadir,subject,filename));
filename = sprintf( '%s_sess-%s_Stim'  ,subject,session); load(fullfile(datadir,subject,filename));
fprintf('plotting.\n')


%%

% GET STIM INFO


% Define TC types
%Freq tuning
%   Freq~=nan
%Level tuning
%   Freq==nan
%   stimDur<800ms
%   numel(unique(dB))>3 ?
%AM rate or depth tuning
%   Freq==nan
%   stimDur>800

% Get all basic tuning stimuli
tStim = Stim(strcmp('basic_tuning',{Stim.stimfile}));
ttrs  = find(strcmp('basic_tuning',{Stim.stimfile}));

Blocks = unique([tStim.block]);

for ib=Blocks
    
    ibStim = tStim([tStim.block]==ib);
    ibTrs  = find([Stim.block]==ib);
    
    if all(~isnan([ibStim.Freq]))                %Freq tuning
%         disp('Freq tuning')
        
        
        
    elseif all(isnan([ibStim.Freq]))...          %Level tuning
            && all([ibStim.stimDur]<800)...
            && numel(unique([ibStim.dB]))>3
%         disp('Level tuning')
        
        
    elseif all(isnan([ibStim.Freq]))...          %AM rate tuning
            && all([ibStim.stimDur]>800)...
            && numel(unique([ibStim.AMrate]))>3
        disp('AM rate tuning')
        
        break
        
        
    elseif all(isnan([ibStim.Freq]))...          %AM depth tuning
            && all([ibStim.stimDur]>800)...
            && numel(unique([ibStim.AMdepth]))>3
%         disp('AM depth tuning')
        
        
    else
        warning('Cannot determine the type of recording from the stimulus parameters.')
        keyboard
        
    end
    
    
end



% Stim struct for the selected tuning curve

Par_matrix = [ [ibStim.AMrate]'...
               [ibStim.dB]'...
               [ibStim.AMdepth]' ] ;
[unique_stim, unique_IDs, StimID] = unique(Par_matrix,'rows','sorted');

% Make raster struct
rateTC = struct(); 
for ks = 1:size(unique_stim,1)
    
    ks_trs = [];
    ks_trs = ibTrs(StimID==StimID(unique_IDs(ks)));
    
    rateTC(ks).block    = ib;
    rateTC(ks).tr_idx   = ks_trs;
    rateTC(ks).AMrate   = unique_stim(ks,1);
    rateTC(ks).dB       = unique_stim(ks,2);
    rateTC(ks).AMdepth  = unique_stim(ks,3);
    rateTC(ks).AMonset  = round(mean([Stim(ks_trs).onsetAM] - [Stim(ks_trs).onset]));
    rateTC(ks).stimDur  = round(mean([Stim(ks_trs).stimDur])) +10; %ramped at end
    rateTC(ks).stim_str = sprintf('AMrate Tuning \n  %i Hz rate, %2.3g dBSPL  \nch %s unit %s',...
                            rateTC(ks).AMrate, rateTC(ks).dB, num2str(channel), num2str(clu) );    
end



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


% Set up raster/histo plot parameters
t_beg  = -499;  %ms
t_end  = 500 + rateTC(ks).stimDur;   %ms
nt     = t_end - t_beg +1;  %each entry 1 ms
bin    = 10;    %ms


smooth.wsize = round(nt/200);   %window size for gaussian smoothing of histo for plotting
smooth.cutoff = 20;   %cutoff for the gaussian smoothing
smooth.stdev = Info.fs/(2*pi*smooth.cutoff); %std for gaussian smoothing of histo for plotting


% Set up figure
nSubPlots = 2;
hS = zeros(numel(rateTC),nSubPlots);

for ks = 1:numel(rateTC)
    
    % Set current figure/subplot handles
    hF(ks) = figure; hold on
    scrsz = get(0,'ScreenSize');
    set(hF(ks),'Position',[1 scrsz(4) scrsz(3) scrsz(4)],...
        'Nextplot','add');
    for isp = 1:nSubPlots
        hS(ks,isp)=subplot(nSubPlots,1,isp);
        set(hS(ks,isp),'Nextplot','add');
    end
    
    rateTC(ks).window_ms = [t_beg t_end];
    
    
    % Get spiketimes for this stim (removing flagged trials)
    tr_this_stim = rateTC(ks).tr_idx;
    
    FLAGGED = [];
    FLAGGED = Info.artifact_trs(channel).trials;
    
    [~,ibt,~] = intersect(tr_this_stim,FLAGGED);
    tr_this_stim(ibt) = [];
    
    
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
    
    rateTC(ks).x = raster_x;
    rateTC(ks).y = raster_y;
    
    hist_raw = sum(hist_raw,1) / it;
    hist_bin = sum(reshape(hist_raw, bin, nt/bin),1)/(bin/1000);
    hist_smooth = smoothts(hist_bin,'g', smooth.wsize, smooth.stdev);
    
    
    
    % Plot this stimulus
    
    suptitle(rateTC(ks).stim_str)
    
    % raster
    subplot(hS(ks,1)); hold on
%     fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*(1+it) + (1+it)/2 ,...
%         Wcolor,'EdgeColor','none')
    plot([0 0],[-30 30],'k:', 'LineWidth', 2)
    plot([rateTC(ks).AMonset rateTC(ks).AMonset],[-30 30],'k:', 'LineWidth', 2)
    plot([rateTC(ks).stimDur rateTC(ks).stimDur],[-30 30],'k:', 'LineWidth', 2)
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
%     fill([1:length(Wave(ks).y) length(Wave(ks).y):-1:1] /Wave(ks).fs*1000, [Wave(ks).y fliplr(-Wave(ks).y)]/4*ymaxval + (ymaxval)/2 ,...
%         Wcolor,'EdgeColor','none')
    plot( t_beg:bin:t_end  , hist_bin , 'k', 'LineWidth', 2)
    plot([0 0],[0 ymaxval],'k:', 'LineWidth', 2)
    plot([rateTC(ks).AMonset rateTC(ks).AMonset],[0 ymaxval],'k:', 'LineWidth', 2)
    plot([rateTC(ks).stimDur rateTC(ks).stimDur],[0 ymaxval],'k:', 'LineWidth', 2)
    set(gca, 'XLim', [t_beg t_end])
    xlabel( 'Time (ms)')
    ylabel('Spikes/sec')
    ylim([0 ymaxval])
    set(gca,'FontSize',figFontSize)
    hold off
    
    set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)
    
    
    
    % SAVE FIGURE
    
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    if ~exist(fullfile(savedir,subject,'^TCs',session),'dir')
        mkdir(fullfile(savedir,subject,'^TCs',session))
    end
    savename = sprintf('%s_%s_ch%i_clu%i_AMrate_%iHz_blk%i',subject,session,channel,clu,rateTC(ks).AMrate,rateTC(ks).block);
%     print(hF(ks),fullfile(savedir,subject,'^TCs',session,savename),'-depsc','-tiff')
    set(gcf,'PaperOrientation','landscape');
    print(hF(ks),'-dpdf',fullfile(savedir,subject,'^TCs',session,savename),'-bestfit')
    
    
    
    % Calculate some basic response properties
    rateTC(ks).nSpk = sum(raster_x > rateTC(ks).AMonset & raster_x < rateTC(ks).stimDur)...
        / it / ((rateTC(ks).stimDur-rateTC(ks).AMonset)/1000); 
        %count spikes during AM, divide by n trials, divide by seconds duration of AM stim
    
    baseRate(ks) = sum(raster_x < 0 )...
        / it / (-t_beg/1000);
    
end


hTC = figure;
plot([rateTC.AMrate],[rateTC.nSpk],':.k','MarkerSize',20)
set(gca,'XScale','log','XTick',[rateTC.AMrate])
hold on
plot([1 256], [mean(baseRate) mean(baseRate)], '--k')
ylim([0 max([rateTC.nSpk])+10])
xlim([1 256])
xlabel('AM rate (Hz)')
ylabel('Firing Rate (spikes/sec)')
unstr = sprintf('%s %s: ch%i clu%i\n',subject,session,channel,clu);
title([unstr 'AM rate Tuning'])


% SAVE TUNING CURVE

savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
savename = sprintf('%s_%s_ch%i_clu%i_TC-AMrate_blk%i',subject,session,channel,clu,rateTC(ks).block);
%     print(hF(ks),fullfile(savedir,subject,'^TCs',session,savename),'-depsc','-tiff')
set(gcf,'PaperOrientation','landscape');
print(hTC,'-dpdf',fullfile(savedir,subject,'^TCs',session,savename),'-bestfit')








end

