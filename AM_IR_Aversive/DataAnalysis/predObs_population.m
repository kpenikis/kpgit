function predObs_population(select_subject, select_session, select_channel, select_clu)
%
%  ap_zFR_population( [subject, session, channel, clu] )
%   All inputs are optional. Any variables not specified will be cycled
%   through.
%   Simply collects the average FR (normalized in z-space) from each trial
%   of periodic stimuli, and plots linear predictions of FR response to IR
%   stimuli against the observed responses.
%   Uses the TrialData (newer) version of saving stimulus info.
%   
%   Excludes datapoints based on: min Ntrials, minFR. Option to exclude MU.
%   Currently allows one unit to be highlighted.
%
%  KP, 2018-03; based on ap_zFR_population from AMstream experiment.
%  updated 2018-04: datapoint filter
%
        
%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  1;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!
AMrespFilt = 1; %add opt to savename/loc
%!!!!!!!!!!!!!!!!!
colorswitch = 'iBMF_FR'; 'iBMF_VS'; 'subject';
%!!!!!!!!!!!!!!!!!
spktimeSHIFT = -0;
%!!!!!!!!!!!!!!!!!
USE_MEASURE = 'FR'; 'TrV'; 'FF';
%!!!!!!!!!!!!!!!!!

switch USE_MEASURE
    case 'FR'
        units = 'z';
        switch units
            case 'hz'
                axmin = -10;
                axmax = 30;
            case 'z'
                axmin = -1;
                axmax = 1.5;
        end
    case {'TrV' 'FF'}
        units = ' ';
        axmin = 0;
        axmax = 4;
end



%% Load Resp struct (with RCorr results)

fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'RespStruct_allSU'));
Resp = q.Resp;


%%
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
vsmallsq = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];
smallsq  = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];

% Prepare the figure
hf1 = figure;
set(hf1,'Position',smallsq,'NextPlot','add')
plot([axmin axmax],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[axmin axmax],'Color',[0.7 0.7 0.7])
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
xlim([axmin axmax])
ylim([axmin axmax])
xlabel('Predicted')
ylabel('Observed')

% Prepare another figure, for individual sequences' observed responses
hf2 = figure;
set(hf2,'Position',smallsq,'NextPlot','add')
plot([axmin axmax],[0 0],'Color',[0.7 0.7 0.7])
hold on
plot([0 0],[axmin axmax],'Color',[0.7 0.7 0.7])
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
xlim([axmin axmax])
ylim([axmin axmax])
xlabel('Predicted')
ylabel('Observed')

subjshapes = {'o' 'o'};

colors = [  0   0   0;...
           84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [colors; 0.7.*bone(2)];
% previous expt : 64 Hz (magenta) : 255 205  60 


%%
% Make empty data table
Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.type     = ' ';
Zdata.UnitN    = 0;
Zdata.IRblock  = 0;
Zdata.spl      = 0;
Zdata.lpn      = 0;
Zdata.zFR_pred = 0;
Zdata.zFR_obs  = 0;

N=0;

%%
% Select a datapoint (or a few?) to highlight

ex_subj = 'XWWf_253400';
ex_sess = 'MA';
ex_ch   = 6;
ex_clu  = 29;


%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    subjects = {select_subject};
else
    subjects = {'WWWf_253400' 'WWWlf_253395'};
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    
    switch subject
        case 'WWWf_253400'
            subjcol = 0.6*[1 1 1];
        case 'WWWlf_253395'
            subjcol = [1 1 1];
    end
     

%%  SESSIONS

% Get list of sessions to check for sorted data

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

% Step through each session
for sess = Sessions'
    
session = char(sess);

if nargin>2 && exist('select_channel','var')
    Channels = select_channel;
else % group data
    
    session = char(sess);
    
    if AMrespFilt %skip sessions with no Robust SUs
        SessResp = Resp(strcmp({Resp.Subject},subject)&strcmp({Resp.Session},session));
        if isempty(SessResp),  continue,   end
        
        % Here add call to a function that outputs the indices of SessResp
        % that contain a responsive unit. With this method, the thresholds can
        % be flexibly changed.
        [sig_ids,SessResp] = identifyResponsiveUnits(SessResp);
        
        % Skip clus that do not show some significant AM response
        Channels = unique([SessResp(sig_ids).Channel]);
        
    else
        Channels =  [1:7 9:16];
    end

end

%%
% Load data files
fn = set_paths_directories(subject,session,1);
fprintf('Loading sess %s...\n',session)
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));

% KS
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],24414)
% clusters = loadKiloSortClusters([fn.sessdata '/sorting'],Info.fs)


%% GET STIM INFO
[dBSPL,LP] = theseSoundParams(TrialData);
AMrates = [2 4 8 16 32];


%% STEP THROUGH EACH CHANNEL

for channel = Channels
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin>3
        Clus = select_clu;
    else
        if AMrespFilt==0
            Clus = select_clu;
            keyboard
        else
            try
                Clus = SessResp([SessResp.Channel]==channel).Clu;
            catch
                continue
            end
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = Clus'
        
        % GET SPIKETIMES
        
        if SUonly && (spikes.labels(spikes.labels(:,1)==clu,2) ~= 2)
            continue
        end
        
        switch spikes.labels(spikes.labels(:,1)==clu,2)
            case 2
                unType = 'SU';
            case 3
                unType = 'MU';
        end
        
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        % spiketrials = spikes.trials(unit_in);
        spiketimes = spiketimes+spktimeSHIFT;
        
        if isempty(spiketimes)
            continue
        end
        
        if nargin<1
        ThisResp = SessResp([SessResp.Channel]==channel & [SessResp.Clu]==clu);
        end
%         fprintf('-----------------------------\n')
        

        
        
        
        
        
        %% ################     NOW THE FUN STARTS     ##################
        
        
        
        
        
        
        
                
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Check the duration of silence at the beginning
        if ((TrialData.offset(1) - TrialData.onset(1))/1000) < 15
            keyboard
        end
        
        bs_hist = 1;
        bs_smth = 20;
        
        [Stream_FRsmooth,Stream_zscore,~,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),bs_hist,bs_smth,'silence');
        
        
        % Check if unit has very low overall firing rate
        if mean(Stream_FRsmooth(TrialData.onset(1):end)) < FRcutoff
            continue
        end
        
        
        
        %%
        
        % Step through each combo of dBSPL, HP, AMdepth
        for spl = dBSPL
            for lpn = LP
                
                fprintf(' analyzing ch %i clu %i\n',channel,clu)
                
                % Get all stimuli presented with these parameters, given a
                % sufficient number of trials without diruptive artifact
                % while the animal was drinking
                
                [all_TDidx,Ntrials,minDur] = get_clean_trials(TrialData,Info.artifact(channel).trials,spl,lpn);
                                
                allStim = unique(TrialData.trID(all_TDidx));
                
                if any(Ntrials(2:6) < minTrs)
                    continue
                end
                
                if sum(Ntrials < minTrs)==1
                    all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<minTrs))  = [];
                    allStim(Ntrials<minTrs)  = [];
                    Ntrials(Ntrials<minTrs) = [];
                elseif  sum(Ntrials < minTrs)>1
                    keyboard
                end
                
                FRmeans = nan(size(allStim'));
                FRvars  = nan(size(allStim'));
                
                for istim = allStim'
                    
                    st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
                    
                    %%%  Skip ITI stimuli 
                    ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
                    
                    for is = 1:numel(ITIflag)
                        TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
                        
                        % Get timestamps of onsets and offsets
                        clear t2 t3 Duration t_win
                        t2 = TrialData.onset(TDidx);
                        t3 = TrialData.offset(TDidx);
                        Duration = mode(diff([t2 t3],1,2));
                        t3 = t2 + Duration;
                        
                        % Now collect FR responses
                        FR_resp = nan(numel(t2),Duration);
                        
                        for it = 1:numel(t2)
                            
                            switch units
                                case 'hz'
                                    FR_resp(it,1:Duration) = Stream_FRsmooth(t2(it):t3(it)-1);
                                case 'z'
                                    FR_resp(it,1:Duration) = Stream_zscore(t2(it):t3(it)-1);
                            end
                            
                        end %it
                        
                        % -- OBSERVED --
                        % Calculate the mean normalized FR for this stimulus
                        FRmeans(istim) = mean(mean(FR_resp,1),'omitnan');
                        FRvars(istim)  = var(mean( FR_resp(:,~isnan(mean(FR_resp,1))) ,2));
                        
                    end %iti filter
                end %istim
                
                
                % -- PREDICTED --
                % Calculate prediction for IR response, based on
                % weighted average of periodic responses
                IR_Prediction = sum( FRmeans(2:6) .* ((1./AMrates)/sum(1./AMrates)) );
                
                % add error bars
                %weighted sum of variances, then converted to overall SEM [--> var(a+b)=var(a)+var(b) ]
                IR_Pred_sem2 = sqrt( sum( FRvars(2:6) .* ((1./AMrates)/sum(1./AMrates)) ) ) / sqrt(length(AMrates));
                %weighted sum of SEMs
                IR_Pred_sem  = sum( ( sqrt(FRvars(2:6))./sqrt(Ntrials(2:6)) ) .* ((1./AMrates)/sum(1./AMrates)) );
                
                
                
                %%
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if nargin>3              % INDIVIDUAL UNIT ZFR MTF
                    
                    xx = [1:5 8+(1:numel(allStim(7:end)))];
                    msize = Ntrials./5;
                    
                    hmtf = figure; hold on
%                     axis square
                    % Plot Warn stim
%                     plot(xx(1),FRmeans(1),'s','MarkerSize',15,...
%                         'MarkerFaceColor','none','MarkerEdgeColor',colors(1,:),'LineWidth',1.5)
%                     plot([xx(1) xx(1)],[FRmeans(1)-(sqrt(FRvars(1))/sqrt(Ntrials(1))) FRmeans(1)+(sqrt(FRvars(1))/sqrt(Ntrials(1)))],...
%                         '-','Color',colors(1,:),'LineWidth',1.5)
                    
                    % Plot periodic stimuli
                    for ir = 1:5
                        plot(xx(ir),FRmeans(ir+1),'o','MarkerSize',15,...
                            'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
                        plot([xx(ir) xx(ir)],[FRmeans(ir+1)-(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1))) FRmeans(ir+1)+(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1)))],...
                            '-','Color',colors(ir+1,:),'LineWidth',2)
                    end
                    
                    % Plot IR prediction
                    plot(7,IR_Prediction,'o','MarkerSize',15,'LineWidth',2,'Color',[32 129 255]./255)
                    plot([7 7],[IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],'-','LineWidth',2,'Color',[32 129 255]./255)
                    
                    % Plot IR observed
                    for ir = 6:length(xx)
                        plot(xx(ir),FRmeans(ir+1),'o','MarkerSize',15,...
                            'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
                        plot([xx(ir) xx(ir)],[FRmeans(ir+1)-(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1))) FRmeans(ir+1)+(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1)))],...
                            '-','Color',colors(ir+1,:),'LineWidth',2)
                    end
                    
                    % Finish formatting
                    set(gca,'XTick',min(xx):max(xx),...
                        'XTickLabel',[Info.stim_ID_key(2:6)' ' ' 'Linear Prediction' ' ' Info.stim_ID_key(allStim(7:end))' ],...
                        'TickLabelInterpreter','none')
                    xtickangle(45)
                    
                    xlim([-1+min(xx) max(xx)+1])
%                     REFERENCE = 'silence';
                    ylabel(sprintf('FR response\n(%s)', units))
                    
                    % Save MTF figure
                    savedir = fullfile(fn.processed,'LinearityFRpopulation');
                    if ~exist(savedir,'dir')
                        mkdir(savedir)
                    end
                    savename = sprintf('%s_%s_ch%i_clu%i_%s_%s_FR-MTF',...
                        subject,session,channel,clu,unType,units);
                    print_eps_kp(hmtf,fullfile(savedir,savename),1)
                    print_svg_kp(hmtf,fullfile(savedir,savename),1)
                    
                    
                    
                    
                    
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                elseif nargin<1          % POPULATION ZFR COMPARISON
                    
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    % Add datapoint to plots
                    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                    
                    N = N+1;
                    
                    % Set color
                    switch colorswitch
                        case 'iBMF_FR'
                            if isfield(ThisResp,'iBMF_FR') && ~isempty(ThisResp.iBMF_FR)
                                plotcol = colors(1+ThisResp.iBMF_FR,:);
                            else
                                plotcol = 'k';
                            end
                        case 'iBMF_VS'
                            if isfield(ThisResp,'iBMF_VS') && ~isempty(ThisResp.iBMF_VS)
                                plotcol = colors(1+ThisResp.iBMF_VS,:);
                            else
                                plotcol = 'k';
                            end
                        case 'subject'
                            plotcol = subjcol;
                    end
                    
                    
                    figure(hf1); hold on
                    
                    % horizontal errorbars (prediction error)
                    plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[mean(FRmeans(7:end)) mean(FRmeans(7:end))],'Color',plotcol)
                    % vertical errorbars (observation error)
                    plot([IR_Prediction IR_Prediction],[mean(FRmeans(7:end))-mean(sqrt(FRvars(7:end))./sqrt(Ntrials(7:end)))...
                        mean(FRmeans(7:end))+mean(sqrt(FRvars(7:end))./sqrt(Ntrials(7:end)))],'Color',plotcol)
                    %                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10))) mean(zFRmeans(7:10))+sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10)))],plotcol)
                    %                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4) mean(zFRmeans(7:10))+sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4)],plotcol)
                    % mean point
                    ip=scatter(IR_Prediction,mean(FRmeans(7:end)),150,subjshapes{subj},'MarkerFaceAlpha',0.45,'MarkerFaceColor',plotcol,'MarkerEdgeColor',plotcol,'LineWidth',1);
                    if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
                        ip.MarkerFaceColor = [32 129 255]./255;
                        ip.MarkerFaceAlpha = 0.65;
                    end
                    
                    figure(hf2); hold on
                    
                    for is = 7:length(allStim)
                        
                        % horizontal errorbars (prediction error)
                        plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[FRmeans(is) FRmeans(is)],'Color',plotcol)
                        % vertical errorbars (observation error)
                        plot([IR_Prediction IR_Prediction],[FRmeans(is)-sqrt(FRvars(is))./sqrt(Ntrials(is)) FRmeans(is)+sqrt(FRvars(is))./sqrt(Ntrials(is))],'Color',plotcol)
                        % mean point
                        ip=scatter(IR_Prediction,FRmeans(is),100,subjshapes{subj},'MarkerFaceAlpha',0.45,'MarkerFaceColor',plotcol,'MarkerEdgeColor',plotcol,'LineWidth',1);
                        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
                            ip.MarkerFaceColor = [32 129 255]./255;
                            ip.MarkerFaceAlpha = 0.65;
                        end
                        
                        
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        % Save data to table
                        % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                        Zdata_addrow = { subject session channel clu unType N is spl lpn IR_Prediction FRmeans(is) };
                        
                        Zdata = [Zdata; Zdata_addrow];
                        
                        
                        
                    end %is
                    
                end %if nargin<1
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
            end %lpn
        end %spl
    end %clu
end %channel
end %sessions
end %subjects

if nargin<1
    
    % Save Z data table
    Zdata(1,:) = [];
    
    
    
    % STATS
    
    p_wsr   = signrank( Zdata.zFR_pred(~isnan(Zdata.zFR_obs)), Zdata.zFR_obs(~isnan(Zdata.zFR_obs)), 'method', 'exact');
    [r,p_r] = corrcoef(Zdata.zFR_pred(~isnan(Zdata.zFR_obs)), Zdata.zFR_obs(~isnan(Zdata.zFR_obs)));
    
    
    % Finish figure (hf2)
    text(axmax/6,2*axmin/4,sprintf('N = %i units\nWSR p=%2.3f\nr=%0.2f, p=%0.2e',N,p_wsr,r(1,2),p_r(1,2)),'FontSize',12)
    
    if SUonly
        if AMrespFilt==0
            title(sprintf('Responses to each IR\nFR (units %s)\nall SUs',units))
            savename = sprintf('%sFR_ObsPred_allSU_eachIR',units);
        else
            title(sprintf('Responses to each IR\nFR (units %s)\nResp SUs',units))
            savename = sprintf('%sFR_ObsPred_RespSU_eachIR',units);
        end
        
    else
        title(sprintf('Responses to each IR\nFR (units %s, ref sil.)\nall Units',units))
        savename = sprintf('%sFR_ObsPred_allUnits_eachIR',units);
    end
    
    % Save population responses figure
    savedir = fullfile(fn.processed,'LinearityFRpopulation');
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    print_eps_kp(hf2,fullfile(savedir,savename))
    print_svg_kp(hf2,fullfile(savedir,savename))
    
    
    
    % Plot histogram of differences between observed and predicted IR
    % responses
    hf1b = figure;
    set(hf1b,'Position',vsmallsq,'NextPlot','add')
    hist(Zdata.zFR_obs - Zdata.zFR_pred,linspace(axmin,axmax,100))
    xlim([axmin -axmin])
    xlabel('Obs-Pred FR, in units of std for each unit')
    title('IR FR responses, Observed - Predicted')
    
    
%     keyboard
    
end %if nargin<1


end %function




