function ap_zFR_population(subject)
%
%  ap_zFR_population(subject)
%    
%
%  KP, 2016-04; last updated 2017-08
%

%  (1,:) = Instantaneous AM rate <-- if Trials stim set, just this
%  (2,:) = Sound output          <-- if Trials stim set, just this
%  (3,:) = AM depth
%  (4,:) = dB SPL
%  (5,:) = HP
%  (6,:) = LP
%  (7,:) = Spout TTL
%  (8,:) = Block label



% IF SAVING PDF FILES
figFontSize      = 14;
rasterMarkerSize = 10;
rasterLineWidth  = 0.5;

% IF SAVING EPS FILES
% figFontSize      = 24;
% rasterMarkerSize = 18;
% rasterLineWidth  = 1;

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)


scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];
figsize2 = [1 scrsz(4) scrsz(3)/2 scrsz(4)];


% Prepare the figure
hf1 = figure;
set(hf1,'Position',figsize1,'NextPlot','add')
plot([-0.5 1],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[-0.5 1],'Color',[0.7 0.7 0.7])
plot([-0.5 1],[-0.5 1],'Color',[0.7 0.7 0.7])
axis square

xlabel('Predicted IR response (zFR)')
ylabel('Observed IR response (zFR)')


% Prepare another figure, for individual sequences' observed responses
hf2 = figure;
set(hf2,'Position',figsize1,'NextPlot','add')
plot([-0.5 1],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[-0.5 1],'Color',[0.7 0.7 0.7])
plot([-0.5 1],[-0.5 1],'Color',[0.7 0.7 0.7])
axis square

xlabel('Predicted IR response (zFR)')
ylabel('Observed IR response (zFR)')

switch subject
    case 'WWWf_244303'
        title('caudal A1')
    case 'WWWr_244300'
        title('DP/VP')
end


histbinsize = 20;
anbinsize   = 50;
smthbinsize = 50;


colors = hsv(6);
colors = [colors; 0.5.*hsv(4)];

unType = {'' 'SU' 'MU'};

N = 0;


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
        
        
        % And set y limits
        ylimvals = [min(Stream_zscore) max(Stream_zscore)];
                
        
        
        %%        
        
        % Also step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                for amd = AMdepth
                    
                    % Get this unit/stimulus combo's N clean blocks (trials) 
                   blocks_N = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    
                    %%
                    % Set up figure
                    fprintf('plotting ch %i clu %i...\n',channel,clu)
%                     hf_t = figure;
%                     set(hf_t,'Position',figsize1,'NextPlot','add')
%                     hf_cv = figure;
%                     set(hf_cv,'Position',figsize1,'NextPlot','add')
                    
                    
                    % Create empty vectors for data
                    FR_pdc_pd = nan(90,4,6);
                    Zresp_pdc = nan(60,6);
                    trtrR = nan(1,10);
                    CV = nan(1,10);
                    for ir = 1:numel(AMrates)
                        eval(sprintf( 'FRresp_IR_%i = [];',AMrates(ir) ))
                        eval(sprintf( 'FRresp_info_%i = [];',AMrates(ir) ))
                    end
                    
                    
                    
                    % Go through each stiulus
                    for ib = 1:numel(Info.blockKey)
                        
                        %if unmodulated or silent, skip for now
                        if ib>=11, continue, end
                        
                        % Get samples of beginning and end of block
                        bkStart_samps = 1+find( diff(SoundData(8,:)==ib) ==  1 );
                        bkStop_samps  =   find( diff(SoundData(8,:)==ib) == -1 );
                        
                        % Remove blocks not this spl, lp noise, am depth,
                        % OR that contain artifact
                        rm_bk = [];
                        for it = 1:numel(bkStart_samps)
                            if ib==5 && (SoundData(8,bkStart_samps(it)-1)==11)
                                rm_bk = [rm_bk it];
                            elseif ib~=5 && (SoundData(8,bkStart_samps(it)-1)==11)
                                keyboard
                            end
                            if ~all(SoundData(4,bkStart_samps(it):bkStop_samps(it))==spl)...               
                                    || ~all(SoundData(6,bkStart_samps(it):bkStop_samps(it))==lpn)...      
                                    || ~all(SoundData(3,bkStart_samps(it):bkStop_samps(it))==amd)...      
                                    || any( intersect(bkStart_samps(it):bkStop_samps(it),ArtifactFlag))...  
                                    || ~all(SoundData(7,bkStart_samps(it):bkStop_samps(it))==1)          
                                rm_bk = [rm_bk it];
                            end
                        end
                        bkStart_samps(rm_bk) = [];
                        bkStop_samps(rm_bk) = [];
                        
                        % Convert to ms
                        bkStart_ms = round( bkStart_samps / Info.fs_sound*1000 );
                        bkStop_ms  = round( bkStop_samps  / Info.fs_sound*1000 );
                        
                        
                        if numel(bkStart_ms) ~= numel(bkStop_ms)
                            keyboard
                        end
                        %view preceeding blocks:  hpb = hist(SoundData(8,bkStart_samps-1),0:12)
                        
                        
                        % Skip if few trials 
                        if blocks_N(ib)<8
                            continue
                        end
                        
                        
                        %% Now prepare to collect response
                        
                        raster_x   = []; 
                        raster_y   = [];
                        stim       = nan(numel(bkStart_ms),max(bkStop_samps-bkStart_samps));
                        hist_raw   = zeros(numel(bkStart_ms),max(bkStop_ms-bkStart_ms));
                        FRresp     = nan(numel(bkStart_ms),1);
                        SpkTs_pdc  = [];
                        SpkTs_pdc2  = [];
                        z_resp     = nan(numel(bkStart_ms),max(bkStop_ms-bkStart_ms));
                        
                        for it = 1:numel(bkStart_ms)
                            
                            % Sound rms envelope
                            stim(it,1:(bkStop_samps(it)-bkStart_samps(it))) = envelope(SoundData(2, bkStart_samps(it):(bkStop_samps(it)-1)),50,'rms');
                            
                            % Get z-score data for this trial
                            z_resp(it,1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_zscore(bkStart_ms(it):bkStop_ms(it)-1);
                            
                            psth_smooth(it,1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_FRsmooth(bkStart_ms(it):bkStop_ms(it)-1);
                            
                            
                            % Get overall spiking data for this block
                            sp=[]; sp = spiketimes( spiketimes>=bkStart_ms(it) & spiketimes<bkStop_ms(it) ) - bkStart_ms(it) +1;
                            hist_raw(it,sp) = 1;
                            raster_x = [raster_x sp];
                            raster_y = [raster_y it .* ones(1,numel(sp))];
                            % Save concatenated spiketimes for VS
%                             SpkTs_pdc2 = [SpkTs_pdc2 sp + (it-1)*(1000/AMrates(ib))];
                            
                            FRresp(it) = numel(sp) / ((bkStop_ms(it)-bkStart_ms(it))/1000);                            
                            
                            % Get spiking data for individual PERIODIC periods
                            if ib<=6 
                                % Estimate samples of period restarts
                                % (phase=0)
                                pd_starts = 0 : round(Info.fs_sound/AMrates(ib)) : (bkStop_samps(it)-bkStart_samps(it));
                                                                
                                % Get spikes
                                which_pd_sec = [0 0.5 1 1.5];
                                for ipd = 1:4
                                    [~,pdc_pd] = min(abs(pd_starts-which_pd_sec(ipd)*Info.fs_sound));
                                    pd_start_ms = round( ( pd_starts(pdc_pd) + bkStart_samps(it) ) /Info.fs_sound *1000 );
                                    FR_pdc_pd(it,ipd,ib) = sum( spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) / (1/AMrates(ib));
                                    SpkTs_pdc = [SpkTs_pdc spiketimes(spiketimes>=pd_start_ms  &  spiketimes<(pd_start_ms + 1000/AMrates(ib)) ) - pd_start_ms + (it-1)*(1000/AMrates(ib))];
                                end
                                
                                % Calculate VS, periodic blocks
                                [VS(ib),RS(ib),P(ib)] = vectorstrength( SpkTs_pdc, AMrates(ib) );
                                
                            end
                            
                            
                            % Get spiking data for each period of IR blocks
                            if ib>6 && ib<11
                                
                                find(cell2mat(cellfun(@eval,strsplit(Info.blockKey{ib}),'UniformOutput',0))==2);
                                
                                newRate = -1+bkStart_samps(it) + find(diff( [SoundData(1, (bkStart_samps(it)-1):(bkStop_samps(it)) ) 0] ));
                                if numel(newRate)~=13
                                    keyboard
                                end
                                for ir = 1:(numel(newRate)-1)
                                    this_rate = SoundData(1, newRate(ir));
                                    t = round( [newRate(ir) newRate(ir+1)-1] / Info.fs_sound*1000 );
                                    sp=[]; sp = spiketimes( spiketimes>=t(1) & spiketimes<=t(2) ) - t(1) +1;
                                    eval(sprintf('FRresp_IR_%i = [FRresp_IR_%i; numel(sp) / (diff(t)/1000)];',this_rate, this_rate ))
                                    eval(sprintf('FRresp_info_%i = [FRresp_info_%i; (newRate(ir)-bkStart_samps(it))/Info.fs_sound ib];', this_rate, this_rate ))
                                end
                                
                                % Calculate VS, IR blocks
                                
                            end
                            
                            
                            
                        end %it
                        
                        
                        
                        % Calculate the normalized mean FR for this block
                        
                        zFRmeans(ib) = mean(mean(z_resp,2,'omitnan'),'omitnan');
                        zFRvars(ib)  = var(mean(z_resp,2,'omitnan'),'omitnan');
                        
                        
                        
                    end %ib
                    
                    
                    % Calculate prediction for IR response, based on
                    % weighted average of periodic responses
                    IR_Prediction = sum( zFRmeans(1:6) .* ((1./AMrates)/sum(1./AMrates)) );
                    
                    % add error bars 
                        % var(a+b) = var(a) + var(b)
                    IR_Pred_sem2  = sqrt(sum(zFRvars(1:6) .* ((1./AMrates)/sum(1./AMrates)))) / sqrt(6);
                    IR_Pred_sem = sum( ( sqrt( zFRvars(1:6))./sqrt(blocks_N(1:6)) ) .* ((1./AMrates)/sum(1./AMrates)) );
                    
                    
                    
                    %% Add datapoint to plots
                    
                    N = N+1;
                    
                    figure(hf1); hold on
                    
                    % horizontal errorbars (prediction error)
                    plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[mean(zFRmeans(7:10)) mean(zFRmeans(7:10))],'k')
                    % vertical errorbars (observation error)
                    plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4) mean(zFRmeans(7:10))+sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4)],'k')
                    % mean point
                    plot(IR_Prediction,mean(zFRmeans(7:10)),'ok','LineWidth',2,'MarkerSize',10,'MarkerFaceColor',[1 1 1])
                    
                    
                    figure(hf2); hold on
                    
                    for ib = 7:10
                        
                    % horizontal errorbars (prediction error)
                    plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[zFRmeans(ib) zFRmeans(ib)],'k')
                    % vertical errorbars (observation error)
                    plot([IR_Prediction IR_Prediction],[zFRmeans(ib)-sqrt(zFRvars(ib))./sqrt(blocks_N(ib)) zFRmeans(ib)+sqrt(zFRvars(ib))./sqrt(blocks_N(ib))],'k')
                    % mean point
                    plot(IR_Prediction,zFRmeans(ib),'ok','LineWidth',2,'MarkerSize',10,'MarkerFaceColor',colors(ib,:))
                    
                    end
                    
                    

                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end % sessions


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




