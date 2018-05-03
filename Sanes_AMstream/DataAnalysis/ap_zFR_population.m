function ap_zFR_population(select_subject,select_session,select_channel,Clus)
%
%  ap_zFR_population(subject)
%    
%
%  KP, 2016-04; last updated 2017-08
%


close all

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!

%%
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];
figsize2 = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

% Prepare the figure
hf1 = figure;
set(hf1,'Position',figsize1,'NextPlot','add')
plot([-1 2],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[-1 2],'Color',[0.7 0.7 0.7])
plot([-1 1.5],[-1 1.5],'Color',[0.7 0.7 0.7])
axis square
xlim([-1 1.5])
ylim([-1 1.5])
xlabel('Predicted')
ylabel('Observed')

% Prepare another figure, for individual sequences' observed responses
hf2 = figure;
set(hf2,'Position',figsize1,'NextPlot','add')
plot([-1 2],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[-1 2],'Color',[0.7 0.7 0.7])
plot([-1 2],[-1 2],'Color',[0.7 0.7 0.7])
axis square
xlim([-1 1.5])
ylim([-1 1.5])
xlabel('Predicted')
ylabel('Observed')

subjshapes = {'o' 'o'};

colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors = [colors; 0.7.*bone(4)];


%%
unType = {'' 'SU' 'MU'};

N = 0;

histbinsize = 20;
anbinsize   = 50;
smthbinsize = 50;


% Make empty data table
Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.UnitN    = 0;
Zdata.IRblock  = 0;
Zdata.stimpars = [0 0];
Zdata.zFR_pred = 0;
Zdata.zFR_obs  = 0;
% Zdata.preUNM   = 0;
% Zdata.postUNM  = 0;


%%
% Set a datapoint (or a few?) to highlight

ex_subj = 'WWWr_244300';
ex_sess = 'IA';
ex_ch   = 1;
ex_clu  = 20;



%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    subjects = {select_subject};
else
    subjects = {'WWWf_244303' 'WWWr_244300'};
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    
    switch subject
        case 'WWWf_244303'
            subjcol = 'k';
        case 'WWWr_244300'
            subjcol = 'k';
    end
     

%%  SESSIONS

% Get list of sessions to check for sorted data

fn = set_paths_directories(subject);

if nargin>1 && exist('select_session','var')
    Sessions = {select_session};
else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
    Sessions = [];
    for ifn = 1:numel(SpkFns)
        if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
            Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
        end
    end
    
end
% Sessions = flipud(Sessions);

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



%% STEP THROUGH EACH CHANNEL

if nargin>2 && exist('select_channel','var')
    Channels = select_channel;
else
    Channels =  [1:7 9:16];
end

for channel = Channels
        
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
            Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = Clus'
        
        % !! Only SU for now !!
        if SUonly && spikes.labels(spikes.labels(:,1)==clu,2) ~= 2
            continue
        end
        
        try
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        
        % Skip units with overall FR in session below a cutoff (3 hz)
        if numel(spiketimes) < round(FRcutoff*length(SoundData)/Info.fs_sound)
            continue
        end
        
                
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Make constant stream of 
        Stream_Spks = zeros(1,1000*ceil((size(SoundData,2)/Info.fs_sound)));
        Stream_Spks(spiketimes) = 1;
        
        %either with standard 20 ms bin
%         Stream_FRbin = 1000*(binspikecounts(Stream_Spks,histbinsize)/histbinsize);
%         Stream_FRbin(isinf(Stream_FRbin)) = nan;
%         foo = repmat(Stream_FRbin,histbinsize,1);
%         Stream_FR = reshape(foo,1,histbinsize*length(Stream_FRbin));
%         Stream_FR = Stream_FR(1:ceil(length(SoundData)/Info.fs_sound*1000));
        
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
        
        REFERENCE = 'unmod';'stream';
        markRED = 0;
        
        switch REFERENCE
            
            case 'stream'
                Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
            
            case 'unmod'
                
                if ~any(strcmp(Zdata.Properties.VariableNames,'preUNM'))
                    Zdata.preUNM   = 0;
                    Zdata.postUNM  = 0;
                end
                
                % Get samples of unmodulated sound and of silence
                unmod_samps  = find(SoundData(8,:)==11);
                
                % Split unmod portion into before and after
                try
                    if ~isempty(find(diff(unmod_samps)>1))
                        unmod1_ms = round( [unmod_samps(1)   unmod_samps(0+find(diff(unmod_samps)>1))] /Info.fs_sound*1000 );
                        unmod2_ms = round( [unmod_samps(1+find(diff(unmod_samps)>1)) unmod_samps(end)] /Info.fs_sound*1000 );
                        
%                         diff(unmod1_ms)/1000
%                         diff(unmod2_ms)/1000
                        
                        meanFR = mean([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        stdFR = std([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        
                    else
                        unmod1_ms = round( [unmod_samps(1) unmod_samps(end)]/Info.fs_sound*1000 );
                        unmod2_ms = [0 0];
                        
                        meanFR = mean(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        stdFR = std(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        
                        fprintf(' !! not much unmod data for sess %s ch%i clu%i \n',session,channel,clu)
%                         markRED = 1;
                    end
                    
                catch
                    keyboard
                end
                
                Stream_zscore = (Stream_FRsmooth(msStart:end) - meanFR) / stdFR;
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
                
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             
        
        
        %%        
        
        % Step through each combo of dBSPL, LP noise, AMdepth
        for spl = dBSPL
            for lpn = LP
                for amd = AMdepth
                    
                    % Get this unit/stimulus combo's N clean blocks (trials)
                    [blocks_N, minTime] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    
                    %%
                    % Set up figure
                    fprintf('plotting ch %i clu %i...\n',channel,clu)
                    
                    
                    % Go through each stiulus
                    for ib = 1:numel(Info.blockKey)
                        
                        %if unmodulated or silent, skip for now
                        if ib>=11, continue, end
                        
                        % Skip if few trials
                        if blocks_N(ib)<minTrs
                            continue
                        end
                        
                        % Get this block start times
                        [~,~,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                            ib,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                        
                        
                        %% Now prepare to collect response
                        
                        z_resp = nan(numel(bkStart_ms),max(bkStop_ms-bkStart_ms));
                        
                        % Collect z-score data for each trial
                        for it = 1:numel(bkStart_ms)
                            
                            z_resp(it,1:(bkStop_ms(it)-bkStart_ms(it))) = Stream_zscore(bkStart_ms(it):bkStop_ms(it)-1);
                             
                        end %it
                        
                        
                        % Calculate the mean normalized FR for this block
                           % now excludes time bins that were not present
                           % for every trial
                        zFRmeans(ib) = mean(mean(z_resp,1),'omitnan');
                        zFRvars(ib)  = var(mean(z_resp(:,~isnan(mean(z_resp,1))),2));
                        
                        
                        
                    end %ib
                    
                    
                    % Calculate prediction for IR response, based on
                    % weighted average of periodic responses
                    IR_Prediction = sum( zFRmeans(1:6) .* ((1./AMrates)/sum(1./AMrates)) );
                    
                    % add error bars 
                        %weighted sum of variances, then converted to overall SEM [--> var(a+b)=var(a)+var(b) ]
                    IR_Pred_sem2 = sqrt( sum( zFRvars(1:6) .* ((1./AMrates)/sum(1./AMrates)) ) ) / sqrt(6);
                        %weighted sum of SEMs
                    IR_Pred_sem  = sum( ( sqrt(zFRvars(1:6))./sqrt(blocks_N(1:6)) ) .* ((1./AMrates)/sum(1./AMrates)) );
                    
                    
                    %%
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    if nargin>3              % INDIVIDUAL UNIT ZFR MTF
                        
                        xx = [1:6 10:13];
                        msize = [((1./AMrates)/sum(1./AMrates)) 0.2 0.2 0.2 0.2];
                        
                        hmtf = figure; hold on
                        axis square
                        for ir = 1:10
                            plot(xx(ir),zFRmeans(ir),'o','MarkerSize',15,...
                                'MarkerFaceColor',colors(ir,:),'MarkerEdgeColor','none')
                            plot([xx(ir) xx(ir)],[zFRmeans(ir)-(sqrt(zFRvars(ir))/sqrt(blocks_N(ir))) zFRmeans(ir)+(sqrt(zFRvars(ir))/sqrt(blocks_N(ir)))],...
                                '-','Color',colors(ir,:),'LineWidth',2)
                        end
                        
                        plot(8,IR_Prediction,'o','MarkerSize',15,'LineWidth',2,'Color',[32 129 255]./255)
                        plot([8 8],[IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],'-','LineWidth',2,'Color',[32 129 255]./255)
                        
                        set(gca,'XTick',1:13,...
                            'XTickLabel',[cellfun(@(x) horzcat(x,' Hz'), Info.blockKey(1:6) ,'UniformOutput',0) ' ' 'Linear Prediction' ' ' Info.blockKey(7:10) ],...
                            'TickLabelInterpreter','none')
                        xtickangle(45)
                        
                        xlim([0 14])
%                         ylim([-0.5 1])
                        ylabel(sprintf('normalized FR response\n(ref''d to %s)', REFERENCE))
                        
                        % Save MTF figure
%                         savedir = fullfile(fn.processed,'zFRpopulation');
                        savedir = fn.anplots;
                        if ~exist(savedir,'dir')
                            mkdir(savedir)
                        end
                        savename = sprintf('%s_%s_ch%i_clu%i_%s_%idB_LP%ihz_ref-%s_zFR-MTF',...
                            subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},spl,lpn,REFERENCE);
                        print_eps_kp(hmtf,fullfile(savedir,savename),1)
                        
                        
                        
                        
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    elseif nargin<1          % POPULATION ZFR COMPARISON
                        
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    % Add datapoint to plots
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    
                    N = N+1;
                    
                    figure(hf1); hold on
                    
                    % horizontal errorbars (prediction error)
                    plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[mean(zFRmeans(7:10)) mean(zFRmeans(7:10))],subjcol)
                    % vertical errorbars (observation error)
                    plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-mean(sqrt(zFRvars(7:10))./sqrt(blocks_N(7:10))) mean(zFRmeans(7:10))+mean(sqrt(zFRvars(7:10))./sqrt(blocks_N(7:10)))],subjcol)
%                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10))) mean(zFRmeans(7:10))+sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10)))],subjcol)
%                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4) mean(zFRmeans(7:10))+sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4)],subjcol)
                    % mean point
                    ip=scatter(IR_Prediction,mean(zFRmeans(7:10)),100,subjshapes{subj},'MarkerFaceAlpha',0.45,'MarkerFaceColor','none','MarkerEdgeColor',subjcol);
                    if markRED
                        ip.MarkerFaceColor = 'r';
                    end
                    if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
                        ip.MarkerFaceColor = [32 129 255]./255;
                        ip.MarkerFaceAlpha = 0.65;
                    end
                    
                    figure(hf2); hold on
                    
                    for ib = 7:10
                        
                        % horizontal errorbars (prediction error)
                        plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[zFRmeans(ib) zFRmeans(ib)],subjcol)
                        % vertical errorbars (observation error)
                        plot([IR_Prediction IR_Prediction],[zFRmeans(ib)-sqrt(zFRvars(ib))./sqrt(blocks_N(ib)) zFRmeans(ib)+sqrt(zFRvars(ib))./sqrt(blocks_N(ib))],subjcol)
                        % mean point
                        ip=scatter(IR_Prediction,zFRmeans(ib),100,subjshapes{subj},'MarkerFaceAlpha',0.45,'MarkerFaceColor','none','MarkerEdgeColor',subjcol);
                        if markRED
                            ip.MarkerFaceColor = 'r';
                        end
                        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
                            ip.MarkerFaceColor = [32 129 255]./255;
                            ip.MarkerFaceAlpha = 0.65;
                        end
                        
                        
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        % Save data to table
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        switch REFERENCE
                            case 'stream'
                                Zdata_addrow = { subject session channel clu N ib [spl lpn] IR_Prediction zFRmeans(ib)};
                                
                            case 'unmod'
                                Zdata_addrow = { subject session channel clu N ib [spl lpn] IR_Prediction zFRmeans(ib) diff(unmod1_ms)/1000 diff(unmod2_ms)/1000};
                        end
                        
                        Zdata = [Zdata; Zdata_addrow];
                        
                        
                        
                    end %ib
                    
                    end %if nargin<1
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    

                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end % sessions

end % subject



%% FINISH AND SAVE FIGURES

if nargin>0
    return
end


% Calculate and print Pearson's r statistic

% For datapoints of each IR block separately 
[r,p]=corrcoef(Zdata.zFR_pred,Zdata.zFR_obs);
figure(hf2)
text(0.75,-0.5,sprintf('r = %4.3f\np = %4.2g',r(2,1),p(2,1)))
text(0.75,-0.2,['N = ' num2str(N)])

% And for datapoints of IR block means
% must recalculate the means from table data
zFR_pred = [];
zFR_obs  = [];
for unN = 1:N
    zFR_pred = [zFR_pred mean(Zdata(Zdata.UnitN==unN,:).zFR_pred)];
    zFR_obs  = [zFR_obs  mean(Zdata(Zdata.UnitN==unN,:).zFR_obs)];
    if (mean(Zdata(Zdata.UnitN==unN,:).zFR_obs) - mean(Zdata(Zdata.UnitN==unN,:).zFR_pred)) > 0.2
        fprintf('obs - pred = %4.3f\n',(mean(Zdata(Zdata.UnitN==unN,:).zFR_obs) - mean(Zdata(Zdata.UnitN==unN,:).zFR_pred)))
        disp(Zdata(Zdata.UnitN==unN,:))
    end
end
[r,p]=corrcoef(zFR_pred,zFR_obs);
figure(hf1)
text(0.75,-0.5,sprintf('r = %4.3f\np = %4.2g',r(2,1),p(2,1)))
text(0.75,-0.2,['N = ' num2str(N)])


% Set save directory

savedir = fullfile(fn.processed,'zFRpopulation');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% Save figs

% mean of IR blocks
savename = sprintf('zFR_ObsPred_SU_%s_IA120',REFERENCE);
set(hf1,'PaperOrientation','landscape');
print_eps_kp(hf1,fullfile(savedir,savename),1)

% each IR sequence plotted out
savename = sprintf('zFR_ObsPred_SU_eachIR_%s_IA120',REFERENCE);
set(hf2,'PaperOrientation','landscape');
print_eps_kp(hf2,fullfile(savedir,savename),1)


%% Finish and save table
% Zdata(1,:)=[];
% 
% savename = sprintf('zFR_ObsPred_SU_%s',REFERENCE);
% writetable(Zdata,fullfile(savedir,savename));



end %function

%}


