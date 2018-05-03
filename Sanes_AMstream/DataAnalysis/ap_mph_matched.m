function ap_mph_matched( select_subject, select_session, channels, clus)
%
%  pp_mph(subject, session, channel, clu)
%    Plots a raster and psth for each unique stimulus. Clu is the label
%    given by UMS (not an index), found in Spikes.sorted.labels.
%

% response metrics: VS, RS, mean phase, nspks, TS?

close all

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10; 
%!!!!!!!!!!!!!!!!!
Classify =  0; 
%!!!!!!!!!!!!!!!!!
Correlations = 1; 
%!!!!!!!!!!!!!!!!!
PlotMPH  =  0;
%!!!!!!!!!!!!!!!!!
plot_which_rates = 1:6;
%!!!!!!!!!!!!!!!!!



% Set some plot options

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',8)

scrsz = get(0,'ScreenSize'); 
figsize1 = [1 scrsz(4) scrsz(3)/3 scrsz(4)];


histbinsize = 20;
anbinsize   = 10;
smthbinsize = 20;

colors1 = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors1 = [colors1; 0.7.*bone(4)];

colors2 = [  0   0  34;...
            50  50  73;...
            94  50  71;...
           159  25  88;...
           222  60  75;...
           244  96  54 ]./255;
colors2 = [colors2; 0.5.*hsv(4)];

colors3 = horzcat(linspace(245,172,6)',linspace(0,245,6)',linspace(66,4,6)')./255;
colors3 = [colors3; 0.7.*bone(4)];

colors4 = [ 51   0   0;...
            50  50  73;...
            94  50  71;...
           159  25  88;...
           222  60  75;...
           244  96  54 ]./255;
colors4 = [colors4; 0.7.*bone(4)];

colors = colors1;

% And create some empty variables

R_IR  = [];
R_pdc = [];

Rdata = table;
Rdata.Subject    = ' ';
Rdata.Session    = ' ';
Rdata.ch         =  0;
Rdata.clu        =  0;
Rdata.AMrate     =  0;
Rdata.stimpars   = [0 0];
Rdata.R_pdc      =  0;
Rdata.R_IR       =  0;
Rdata.diff_pdcIR =  0;
Rdata.R_context  = {magic(3)};





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
    
    % Set paths
    fn = set_paths_directories(subject);
    savedir = fullfile(fn.processed,subject,'^an_plots','MPH');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end


%%  SESSIONS

if ~exist('select_session','var')
    
% Get list of sessions to check for sorted data

SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
Sessions = [];
for ifn = 1:numel(SpkFns)
    if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
        Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
    end
end

else
    Sessions = {select_session};
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
%             disp(' SESSION MAY NOT BE MANUALLY SORTED YET')
%             return
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
                    % Create empty Table for data
                    
                    MPH = table;
                    
                    MPH.AMrate    = 0;
                    MPH.Context   = 'placeholder'; 
                    MPH.SeqPos    = nan;
                    MPH.Starttime = nan;
%                     MPH.PrevBlock = nan;
%                     MPH.PrevRate  = nan;
                    MPH.raster    = {magic(3)};
                    MPH.stimpars  = [0 0];
                    
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % First collect data from IR periods
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for ib = 7:10
                        
                        % Create empty struct for gathering data
                        MPH_temp = struct();
                        for ir = 1:numel(AMrates)
                            MPH_temp(ir).raster  = zeros(2*min(blocks_N),ceil(1000/AMrates(ir)));
                        end
                        
                        % Get onsets for this block type
                        [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                            ib,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                        
                        % Permute trials to pick random number
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
                                
                                % Determine whether this period was in the first or second sequence of the IR block
                                seqPos = 1+ceil( (newRates(ir)-bkStart_samps(trials(it))) /Info.fs_sound - sum(1./AMrates)+0.005 );
                                if (seqPos<1 || seqPos>2)
                                    keyboard
                                end
                                
                                % Save some info about this period
                                MPH_temp(AMrates==this_rate).pdtime(min(blocks_N)*(seqPos-1)+it,1)  = seqPos;
                                MPH_temp(AMrates==this_rate).pdtime(min(blocks_N)*(seqPos-1)+it,2)  = 1000*(newRates(ir)-bkStart_samps(trials(it))) /Info.fs_sound;
                                MPH_temp(AMrates==this_rate).pdtime(min(blocks_N)*(seqPos-1)+it,3)  = SoundData(8, bkStart_samps(trials(it))-1 );
                                MPH_temp(AMrates==this_rate).pdtime(min(blocks_N)*(seqPos-1)+it,4)  = SoundData(1, bkStart_samps(trials(it))-1 );
                                
                                % Put spikes into corresponding raster
                                MPH_temp(AMrates==this_rate).raster(min(blocks_N)*(seqPos-1)+it,sp) = 1; 
                                
                                
                            end
                            
                        end %it
                        
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        % Put the data just collected into the table
                        %    [ AMrate Context SeqPos Starttime PrevBlock PrevRate raster ]
                        for ir = 1:numel(AMrates)
                            for isp = unique(round(MPH_temp(ir).pdtime(:,1)))'
                                
                                these_rows = find(round(MPH_temp(ir).pdtime(:,1))==isp);
                                
                                context   = ['IR_' num2str(ib)];
                                starttime = mode(round(MPH_temp(ir).pdtime(these_rows,2)));
                                prevblock = mode(round(MPH_temp(ir).pdtime(these_rows,3)));
                                prevrate  = mode(round(MPH_temp(ir).pdtime(these_rows,4)));
                                raster    = MPH_temp(ir).raster(these_rows,1:ceil(1000/AMrates(ir)));
                                
                                MPH_addrow = {AMrates(ir) context isp starttime {raster} [spl lpn]};
                                MPH = [MPH; MPH_addrow];
                                
                            end
                        end
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                    end %ib
                    clear MPH_temp
                    MPH(1,:) = [];
                    
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Now collect data from periodic blocks
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    for irate = 1:6
                        
                        % Get period starttimes to match
                        
                        starttimes_to_match = sort(MPH(MPH.AMrate==AMrates(irate),:).Starttime);
                        pd_starts = 0 : 1000/AMrates(irate) : 1999 ;
                        
                        which_pds = [];
                        for ii=1:size(starttimes_to_match,1)
                            [~,idx] = min(abs(pd_starts-starttimes_to_match(ii)));
                            which_pds = [which_pds pd_starts(idx)];
                        end
                        
                        which_pds = unique(which_pds)/1000;
                        
                        % Create empty struct for gathering data
                        MPH_temp = struct();
                        for ipd = 1:numel(which_pds)
                            MPH_temp(ipd).raster  = zeros(min(blocks_N),ceil(1000/AMrates(irate)));
                        end
                        
                        
                        % Get onsets for this block type
                        [bkStart_samps,bkStop_samps,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                            irate,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                        
                        trials = randperm(numel(bkStart_ms));
                        
                        for it = 1:numel(trials)
                            if it>min(blocks_N), continue, end
                            
                            % Estimate samples of period restarts
                            % (phase=0)
                            pd_starts = 0 : round(Info.fs_sound/AMrates(irate)) : (bkStop_samps(trials(it))-bkStart_samps(trials(it)));
                            
                            for ipd = 1:length(which_pds)
                                
                                [~,pdc_pd] = min(abs(pd_starts-which_pds(ipd)*Info.fs_sound));
                                pd_start_ms = round( ( pd_starts(pdc_pd) + bkStart_samps(trials(it)) ) /Info.fs_sound *1000 );
                                
                                MPH_temp(ipd).pdtime(it,1)  = ipd;
                                MPH_temp(ipd).pdtime(it,2)  = round( pd_starts(pdc_pd)/Info.fs_sound *1000 );
                                
                                % Get spikes and put into corresponding raster
                                sp=[]; sp = spiketimes(  spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(irate)) ) - pd_start_ms +1 ;
                                MPH_temp(ipd).raster(it,sp) = 1;
                                
                            end
                            
                        end %it
                        
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        % Put the data just collected into the table
                        %    [ AMrate Context SeqPos Starttime raster ]
                        for ipd = 1:numel(which_pds)
                                
                                context   = 'periodic';
                                starttime = mode(round(MPH_temp(ipd).pdtime(:,2)));
                                raster    = MPH_temp(ipd).raster(:,1:ceil(1000/AMrates(irate)));
                                
                                MPH_addrow = {AMrates(irate) context ipd starttime {raster} [spl lpn]};
                                MPH = [MPH; MPH_addrow];
                                
                        end
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                    end %irate
                    
                    
                    
                    
                    
                    
                    %% MPH plots
                                        
                    for ir = plot_which_rates
                        
                        
                        Data_IR  = MPH( MPH.AMrate==AMrates(ir) & strncmp(MPH.Context,'IR',2) ,:);
                        Data_pdc = MPH( MPH.AMrate==AMrates(ir) & strncmp(MPH.Context,'per',3) ,:);
                        
                        [~,row_idx] = sort(Data_IR.Starttime);
                        
                        
                        %% Compare MPH correlations
                        
                        if Correlations
                            
                            % Comparing COLUMNS
                            % ((Are temporal responses in one context more consistent than the other?))
                            [meanR_IR,meanR_pdc] = corr_comp_mph(Data_IR,Data_pdc);
                            
                            R_IR  = [R_IR;  meanR_IR];
                            R_pdc = [R_pdc; meanR_pdc];
                            
                            
                            % Comparing ROWS
                            % ((Is the temporal response similar across contexts?))
                            R_context = corr_context_mph(Data_IR,Data_pdc,AMrates(ir));
                            
                            
                            % Add data to table
                            Rdata_addrow = {subject session channel clu ir Data_IR(1,:).stimpars meanR_pdc meanR_IR meanR_pdc-meanR_IR {R_context}};
                            Rdata = [Rdata; Rdata_addrow];
                                                        
                        end
                        
                        
                        %%
                        if PlotMPH
                            
                            % Make figure
                            hf(ir) = figure;
                            set(hf(ir),'Position',figsize1,'NextPlot','add')
                            
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
                            
                            for irow = 1:numel(row_idx)
                                
                                sp_idx = 2*irow + 2;
                                
                                %~~~~~~~~~~
                                % IR first
                                
                                % Get raster data
                                raster_x = []; raster_y = [];
                                for it = 1:size(Data_IR(row_idx(irow),:).raster{:},1)
                                    raster_x = [raster_x find(Data_IR(row_idx(irow),:).raster{:}(it,:))];
                                    raster_y = [raster_y repmat(nrows+it,1,numel(find(Data_IR(row_idx(irow),:).raster{:}(it,:))))];
                                end
                                
                                % Plot rasters (left early, right middle)
                                ax_r(2)=subplot(9,2,2);
                                hold on; box off
                                hL=plot(raster_x, raster_y,'+','Color',contextcolors(1+find(strcmp(Data_IR(row_idx(irow),:).Context,unique(Data_IR.Context,'stable'))),:),...
                                    'MarkerSize',2,'LineWidth',0.5);
                                title(num2str(meanR_IR))
                                
                                % Plot folded MPH
                                ax_m(sp_idx-2) = subplot(9,2,sp_idx);
                                sphist = reshape([hist(raster_x,linspace(0,1000/AMrates(ir),52)); hist(raster_x,linspace(0,1000/AMrates(ir),52))],1,52*2);
                                patch(reshape([0.5:52.5; 0.5:52.5],1,53*2), [0 sphist 0], contextcolors(1+find(strcmp(Data_IR(row_idx(irow),:).Context,unique(Data_IR.Context,'stable'))),:) )
                                box off
                                title([Data_IR(row_idx(irow),:).Context{:} ', ' num2str(Data_IR(row_idx(irow),:).Starttime) 'ms'])
                                
                                % Keep track of largest bin spk count to
                                % set ymax
                                if max(hist(raster_x,linspace(0,1000/AMrates(ir),52))) > max_count
                                    max_count = max(hist(raster_x,linspace(0,1000/AMrates(ir),52)));
                                end
                                
                                
                                %~~~~~~~~~~~~~~
                                % Now Periodic
                                
                                [~,this_pdc] = min(abs(Data_pdc.Starttime - Data_IR(row_idx(irow),:).Starttime));
                                
                                % Get raster data
                                raster_x = []; raster_y = [];
                                for it = 1:size(Data_pdc(this_pdc,:).raster{:},1)
                                    raster_x = [raster_x find(Data_pdc(this_pdc,:).raster{:}(it,:))];
                                    raster_y = [raster_y repmat(nrows+it,1,numel(find(Data_pdc(this_pdc,:).raster{:}(it,:))))];
                                end
                                
                                % Plot rasters (left early, right middle)
                                ax_r(1)=subplot(9,2,1);
                                hold on; box off
                                hL=plot(raster_x, raster_y,'+','Color',contextcolors(1,:),...
                                    'MarkerSize',2,'LineWidth',0.5);
                                xlabel('time (ms)')
                                ylabel([num2str(min(blocks_N)) ' trs each'])
                                title(num2str(meanR_pdc))
                                
                                % Plot folded MPH
                                ax_m(sp_idx-2-1) = subplot(9,2,sp_idx-1); box off
                                sphist = reshape([hist(raster_x,linspace(0,1000/AMrates(ir),52)); hist(raster_x,linspace(0,1000/AMrates(ir),52))],1,52*2);
                                patch(reshape([0.5:52.5; 0.5:52.5],1,53*2), [0 sphist 0], contextcolors(1,:) )
                                box off
                                title([Data_pdc(this_pdc,:).Context{:} ', ' num2str(Data_pdc(this_pdc,:).Starttime) 'ms'])
                                
                                % Keep track of largest bin spk count to
                                % set ymax
                                if max(hist(raster_x,linspace(0,1000/AMrates(ir),52))) > max_count
                                    max_count = max(hist(raster_x,linspace(0,1000/AMrates(ir),52)));
                                end
                                
                                
                                
                                % Keep track of rows for raster plot
                                nrows = nrows+min(blocks_N);
                                
                                
                                
                            end %irow
                            
                            
                            
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
                            
                            
                            % SAVE MPH FIGURE
                            
                            savename = sprintf('%s_%s_ch%i_clu%i_%s_%idB_LP%ihz_MPH_%iHz',...
                                subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn,AMrates(ir));
                            set(hf(ir),'PaperOrientation','landscape');
                            print(hf(ir),'-dpdf',fullfile(savedir,savename),'-bestfit')
                            print(hf(ir),fullfile(savedir,savename),'-depsc','-tiff')
                            
                            
                        end %if PlotMPH
                        
                        
                        
                    end %ir
                    
                    
                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
        
end %channel

end %session

end %subject

if nargin>1
    return
end


%% Finish and save table

Rdata(1,:) = [];

% Save
% savename = sprintf('MPH-corr_pdcIR_%s',subject);
savedir = fullfile(fn.processed,'MPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
writetable(Rdata,fullfile(savedir,'MPH-corr_pdcIR'));




%%  MPH shape similarity within context

% Stats

% Get correlation of average correlation coefficients
[r,p_r]=corrcoef(Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_pdc, Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_IR);
% Perform t-test of periodic/IR correlation coefficients
[t,p_t]=ttest2(Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_pdc, Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_IR);
% Perform non-parametric Wilcoxon-Mann_Whitney test of periodic/IR correlation coefficients
p_wsr=signrank(Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_pdc, Rdata(~isnan(Rdata.R_pdc)&~isnan(Rdata.R_IR),:).R_IR);


% Plot 
set(0,'DefaultAxesFontSize',14)

hfr=figure; 
figsize2 = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];
set(gcf,'Position',figsize2)
plot([-0.1 1],[-0.1 1],'k-')
hold on
axis square
xlim([-0.1 1])
ylim([-0.1 1])
% title(sprintf('pairwise correlations of \nMPHs constructed from \nperiodic vs IR contexts'))
xlabel('avg R - Periodic')
ylabel('avg R - IR')
set(gca,'xtick',[],'ytick',[])

subjshapes = {'o' 'o'};
for subj = 1:numel(subjects)
    subject = subjects{subj};
    
    for ir = unique(Rdata.AMrate)'
        
        scatter( Rdata((Rdata.AMrate==ir) & (strcmp(Rdata.Subject,subject)),:).R_pdc,...
            Rdata((Rdata.AMrate==ir) & (strcmp(Rdata.Subject,subject)),:).R_IR,...
            60,subjshapes{subj}, 'MarkerEdgeColor','k',...
            'MarkerFaceAlpha',0.65,...
            'MarkerFaceColor','none' )
         
    end
end
text(0.8,0.1,sprintf('r = %4.3f\np = %4.2g',r(2,1),p_r(2,1)))
text(0.7,0.1,sprintf('p = %4.3f',p_wsr))


% Save figure
savename = 'MPH_ContextSimilarityComp_bw';
print_eps_kp(hfr,fullfile(savedir,savename),1);


%%  MPH shape similarity across contexts

% 
% % Collect data
% R_context = nan(size(Rdata,1),8);
% for idp = 1:size(Rdata,1)
%     R_context(idp,:) = Rdata(idp,:).R_context{:}(:,1)';
% end
% 
% % Stats
% 
% % Plot 
% hfr=figure; 
% figsize2 = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
% set(gcf,'Position',figsize2)
% for ipd = 1:8
%     hsp(ipd)=subplot(8,1,ipd);
%     hist(R_context(:,ipd),-0.2:0.05:1)
% end
% linkaxes(hsp,'xy')
% ylim([0 30])
% xlim([-0.3 1])
% suptitle(sprintf('Distribution of correlations between IR MPH\nand corresponding PDC period,\neach IR period separated'))
% 
% 
% R_context_vector = reshape(R_context,1,size(R_context,1)*size(R_context,2));
% R_context_vector(~isnan(R_context_vector));
% 
% hfr=figure; 
% figsize2 = [1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4];
% set(gcf,'Position',figsize2)
% hist(R_context_vector(~isnan(R_context_vector)),-0.2:0.05:1)
% title('Distribution of correlations between IR MPH and its corresponding PDC period, all data aggregated')
% xlim([-0.3 1])
% 
% 
% %  now separating by AM rate instead of period #
% 
% % Collect data
% R_context = nan(size(Rdata,1),8,6);
% for idp = 1:size(Rdata,1)
%     R_context(idp,:,Rdata(idp,:).AMrate) = Rdata(idp,:).R_context{:}(:,1)';
% end
% R_context_mat = reshape(R_context,[size(R_context,1)*size(R_context,2),6]);
% 
% hfr=figure; 
% figsize2 = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
% set(gcf,'Position',figsize2)
% for ir = 1:6
%     hsp(ir)=subplot(6,1,ir);
%     hist(R_context_mat(:,ir),-0.2:0.05:1)
% end
% linkaxes(hsp,'xy')
% xlim([-0.3 1])
% suptitle(sprintf('Distribution of correlations between IR MPH\nand corresponding PDC period,\neach AM rate separated'))
% 
% 
% 


%%

% Rates for which the periodic MPH is more consistent than IR
Rdata( (Rdata.R_pdc>Rdata.R_IR),:);

% Rates for which the IR MPH is more consistent than periodic
Rdata( (Rdata.R_pdc<Rdata.R_IR),:);



end %function


