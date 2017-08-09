function Latency = ap_plot_TC_AMrate(subject, session, channel)%, clu, HP_, LP_, dB_ )
%
%  ap_plot_TC_AMrate(subject, session, channel)  
%    Tuning analyses for AM stream experiment (freq and level/noiseband).
%    Plots a raster and psth for each unique stimulus. Then plots TC for
%    the entire block.
%    Also calculates latency if the FR goes significantly above baseline.
%    
%
%  KP, 2017-07


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
close all

Latency = [];

%%%
ymaxval = 150;
%%%


% Load data files
fn = set_paths_directories(subject,session);
fprintf('loading data...\n')
filename = sprintf( '%s_sess-%s_Spikes',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'  ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Stim'  ,subject,session); load(fullfile(fn.processed,subject,filename));

fprintf('plotting.\n')


%%
% GET SPIKE TIMES
spikes = Spikes.sorted(channel);

% Find clus to plot
if nargin<4
    if all(spikes.labels(:,2)==1)
        disp(' SESSION MAY NOT BE MANUALLY SORTED YET')
        return
    end
    if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
        disp(' !! no valid clus for this channel')
        return
    else
        clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3,1);
    end
end

% STEP THROUGH EACH CLU

for clu=clus'
    
unit_in = find(spikes.assigns==clu);
spiketimes = round(spikes.spiketimes(unit_in) * 1000);  %ms
spiketrials = spikes.trials(unit_in); 

if isempty(spiketimes)
    error('no spike events found for this clu')
elseif spikes.labels(spikes.labels(:,1)==clu,2) == 4
    warning('  this clu is labeled as noise. are you sure you want to plot?')
    keyboard
end


%%
% GET STIM INFO

Blocks = unique([Stim.block]);
Resp=[];

% Go through each recording block
for ib=Blocks
    
    ibStim = Stim([Stim.block]==ib);
    ibTrs  = find([Stim.block]==ib);
    
    % First determine what type of tuning was preseted this block
    
    if all(~isnan([ibStim.Freq])) ...
            && all(isnan([ibStim.LP]))           %Freq tuning
        disp('Freq tuning')
        
        TCtype = 'Freq';
        Par_matrix = [ [ibStim.Freq]'...
               [ibStim.dB]' ];
        
        
    elseif all(isnan([ibStim.Freq]))...          %Level tuning
            && numel(unique([ibStim.dB]))>3
        disp('Level tuning')
        
        TCtype = 'Noise';
        Par_matrix = [ [ibStim.LP]'...
               [ibStim.dB]'...
               [ibStim.HP]' ] ;
        
    else
        warning('Cannot determine the type of recording from the stimulus parameters.')
        keyboard
        
    end
    
    
    
    % Make struct for the selected tuning curve
    
    [unique_stim, unique_IDs, StimID] = unique(Par_matrix,'rows','sorted');
    
    stimTC = struct();
    for ks = 1:size(unique_stim,1)
        
        ks_trs = [];
        ks_trs = ibTrs(StimID==StimID(unique_IDs(ks)));
        ks_trs([Stim(ks_trs).stimDur]<0) = [];
        
        stimTC(ks).block    = ib;
        stimTC(ks).TCtype   = TCtype;
        stimTC(ks).tr_idx   = ks_trs;
        stimTC(ks).dB       = unique_stim(ks,2);
        stimTC(ks).stimDur  = round(mean([Stim(ks_trs).stimDur])) +10; %ramped at end
        
        if strcmp(TCtype,'Freq')
            stimTC(ks).Freq     = unique_stim(ks,1);
            stimTC(ks).HP       = nan;
            stimTC(ks).LP       = nan;
            stimTC(ks).stim_str = sprintf('%s Tuning\n%i Hz, %i dB',TCtype,stimTC(ks).Freq, stimTC(ks).dB);
        else
            stimTC(ks).LP       = unique_stim(ks,1);
            stimTC(ks).HP       = unique_stim(ks,3);
            stimTC(ks).Freq     = nan;
            stimTC(ks).stim_str = sprintf('%s Tuning\n%i-%i Hz, %i dB',TCtype,stimTC(ks).HP, stimTC(ks).LP, stimTC(ks).dB);
        end
        
    end



%%

% Set up raster/histo plot parameters
t_beg  = -199;  %ms
t_end  = 200 + stimTC(ks).stimDur;   %ms
nt     = t_end - t_beg +1;  %each entry 1 ms
bin    = 10;    %ms


smooth.wsize = round(nt/200);   %window size for gaussian smoothing of histo for plotting
smooth.cutoff = 20;   %cutoff for the gaussian smoothing
smooth.stdev = Info.fs/(2*pi*smooth.cutoff); %std for gaussian smoothing of histo for plotting


% Set up figure
nSubPlots = 2;
hS = zeros(numel(stimTC),nSubPlots);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% STEP THROUGH EACH STIMULUS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ks = 1:numel(stimTC)
    
    % Set current figure/subplot handles
    hF(ks) = figure; hold on
    scrsz = get(0,'ScreenSize');
    set(hF(ks),'Position',[1 scrsz(4) scrsz(3) scrsz(4)],...
        'Nextplot','add');
    for isp = 1:nSubPlots
        hS(ks,isp)=subplot(nSubPlots,1,isp);
        set(hS(ks,isp),'Nextplot','add');
    end
    
    stimTC(ks).window_ms = [t_beg t_end];
    
    
    % Get trials for this stim, and remove flagged ones
    tr_this_stim = stimTC(ks).tr_idx;
    
    FLAGGED = [];
    FLAGGED = Info.artifact_trs(channel).trials;
    
    [~,ibt,~] = intersect(tr_this_stim,FLAGGED);
    tr_this_stim(ibt) = [];
    
    % Skip stimuli with fewer than 8 clean trials
    if numel(tr_this_stim)<8
        continue
    end
    
    % Get spiking data
    raster_x=[];  raster_y=[];  sp_raw=zeros(1,nt);
    for it = 1:numel(tr_this_stim)
        
        sp=[];  spk_in=[];
        spk_in = find(spiketrials==tr_this_stim(it));
        sp = spiketimes(spk_in) + ones(size(spiketimes(spk_in)))*(Info.t_win_ms(1)-1); %ms, rel to t0
        sp = sp( sp>=t_beg & sp<= t_end );
        
        sp_raw(it, sp-t_beg+1) = 1;
        raster_x = [raster_x sp];
        raster_y = [raster_y it .* ones(1,numel(sp))];
        nsp(it)  = sum(sp<0)/-t_beg*1000;
    end
    
    stimTC(ks).x = raster_x;
    stimTC(ks).y = raster_y;
    
    hist_raw = sum(sp_raw,1) / it;
    hist_bin = sum(reshape(hist_raw, bin, nt/bin),1)/(bin/1000);
%     hist_smooth = smoothts(hist_bin,'g', smooth.wsize, smooth.stdev);
    hist_smooth = smoothFR(sp_raw);
    
    
    
    % Calculate some basic response properties
    
    stimTC(ks).nSpk = sum(raster_x > 0 & raster_x < stimTC(ks).stimDur)...
        / it / (stimTC(ks).stimDur/1000); 
        %count spikes during stim, divide by n trials, divide by seconds duration of stim
    
    stimTC(ks).baseRate = sum(raster_x < 0 )...
        / it / (-t_beg/1000);
    
    
    %%%%%
    
    % Calculate latency if FR went significantly above baseline
    
    thresh = mean(nsp)+2*std(nsp);
    
    stimTC(ks).latency = nan;
    if any(hist_smooth >thresh )
        signf = find(hist_smooth>thresh);
        signf(signf<-t_beg)=[];
        if ~isempty(signf)
            stimTC(ks).latency = signf(1)+t_beg;
        end
    end
    
    % Save latency as output variable if this is the matching stimulus
%     if ( strcmp(TCtype,'Noise')...
%             && stimTC(ks).HP == HP_ ...
%             && stimTC(ks).LP == LP_ ...
%             && stimTC(ks).dB == dB_ )
%         
%         keyboard
%         Latency = stimTC(ks).latency;
%     end
    
    %%%%%%
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Plot this stimulus
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    suptitle(stimTC(ks).stim_str)
    
    % raster
    subplot(hS(ks,1)); hold on
    plot([0 0],[-30 30],'k:', 'LineWidth', 2)
    plot([stimTC(ks).stimDur stimTC(ks).stimDur],[-30 30],'k:', 'LineWidth', 2)
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
%     plot( t_beg:bin:t_end  , hist_bin , 'k', 'LineWidth', 2)
    plot( t_beg:t_end, hist_smooth , 'k', 'LineWidth', 2)
    plot([0 0],[0 ymaxval],'k:', 'LineWidth', 2)
    plot([stimTC(ks).stimDur stimTC(ks).stimDur],[0 ymaxval],'k:', 'LineWidth', 2)
    plot(stimTC(ks).latency, hist_smooth(max(stimTC(ks).latency-t_beg,1)),'g*','MarkerSize',20,'LineWidth', 2)
    set(gca, 'XLim', [t_beg t_end])
    xlabel( 'Time (ms)')
    ylabel('Spikes/sec')
    ylim([0 ymaxval])
    set(gca,'FontSize',figFontSize)
    hold off
    
    set(findall(gcf, 'Type','text'), 'FontSize', figFontSize)
    
    
    
    % SAVE FIGURE
    
%     savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
%     if ~exist(fullfile(savedir,subject,'^TCs',session),'dir')
%         mkdir(fullfile(savedir,subject,'^TCs',session))
%     end
%     savename = sprintf('%s_%s_ch%i_clu%i_AMrate_%iHz_blk%i',subject,session,channel,clu,stimTC(ks).AMrate,stimTC(ks).block);
% %     print(hF(ks),fullfile(savedir,subject,'^TCs',session,savename),'-depsc','-tiff')
%     set(gcf,'PaperOrientation','landscape');
%     print(hF(ks),'-dpdf',fullfile(savedir,subject,'^TCs',session,savename),'-bestfit')
%     
    
    
end


%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% MAKE TUNING CURVE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hTC = figure;
clear xdata ydata legstr

% Prepare data
if strcmp(TCtype,'Freq')
    if numel(unique([stimTC.dB])) >1
        keyboard
    end
    xdata = [stimTC.Freq];
    ydata = [stimTC.nSpk];
    xlabelstr = 'Frequency (Hz)';
    xscaletype = 'log';
    legstr = [num2str(stimTC(1).dB) 'dB'];

else
    xlabelstr = 'dB SPL';
    xscaletype = 'linear';
        for np = unique([stimTC.LP])
            xdata(np==unique([stimTC.LP]),:) = [stimTC([stimTC.LP]==np).dB];
            ydata(np==unique([stimTC.LP]),:) = [stimTC([stimTC.LP]==np).nSpk];
            legstr{np==unique([stimTC.LP])} = num2str(np);
        end
end

% Plot it
for il=1:size(xdata,1)
    hold on
    hp(il)=plot(xdata(il,:),ydata(il,:),':.','MarkerSize',20);
end
set(gca,'XScale',xscaletype,'XTick',xdata(1,:))
legend(hp,legstr)
plot([1 max(xdata(1,:))+max(xdata(1,:))/2], [mean([stimTC.baseRate]) mean([stimTC.baseRate])], '--k')
ylim([0 max([stimTC.nSpk])+10])
xlim([min(xdata(1,:))-min(xdata(1,:))/2 max(xdata(1,:))+max(xdata(1,:))/2])
xlabel(xlabelstr)
ylabel('Firing Rate (spikes/sec)')
unstr = sprintf('%s %s: ch%i clu%i\n',subject,session,channel,clu);
title([unstr TCtype])


% SAVE TUNING CURVE
if ~exist(fn.tuning,'dir')
    mkdir(fn.tuning)
end
savename = sprintf('%s_%s_ch%i_clu%i_blk%i_TC-%s',subject,session,channel,clu,stimTC(ks).block,TCtype);
set(gcf,'PaperOrientation','landscape');
print(hTC,'-dpdf',fullfile(fn.tuning,savename),'-bestfit')



Resp = [Resp stimTC];

end % ib=Blocks


% Save Resp struct 
% fprintf('\nsaving data...\n')
% savename = sprintf('%s_sess-%s_ch%i_clu%i_Resp',subject,session,channel,clu);
% save(fullfile(fn.processed,subject,savename),'Resp','-v7.3');


end % clu=clus

end

