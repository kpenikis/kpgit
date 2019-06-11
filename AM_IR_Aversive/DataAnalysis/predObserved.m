function predObserved(USE_MEASURE)
% predObserved(USE_MEASURE)
% 
%   Based on MTF of each unit for Pdc stim, calculate predicted IR response.
%   Compare this prediction to the observed response, for each unit in the
%   Units files. Highlight the ones that are statistically different than
%   prediction. Plot MTF of an example unit if desired.
% 
% KP, 2018-07-01
%  based on predObs_population and previous code, but now uses Unit files
%


% for significance calculation must measure each IR stimulus separately
% (add in fig 2 and run again)


close all
global fn AMrates rateVec_AC rateVec_DB 

%!!!!!!!!!!!!!!!!!
exclOnset   = 0; 
%!!!!!!!!!!!!!!!!!
plotMTF     = 1;
%!!!!!!!!!!!!!!!!!
minTrs      = 10;
%!!!!!!!!!!!!!!!!!
colorswitch = 'subject'; 'iBMF_FR'; 'iBMF_VS'; 
%!!!!!!!!!!!!!!!!!
alfa        = 0.05;
%!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


if nargin<1
    USE_MEASURE =  'FF'; 'FR'; 'TrV';
end

switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 10^-2;
                axmax = 10^2; %10*ceil(2+max([UnitData.BaseFR])/10);
                axscale = 'log';
            case 'z'
                axmin = -1;
                axmax = 1.5;
                axscale = 'linear';
        end
    case {'TrV' 'FF'}
        units = 'Hz';
        axmin = 0.01;
        axmax = 6;
        axscale = 'linear';
end



%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');
vsmallsq = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];
smallsq  = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];


%~~~~~~~~~~~~~~~~~~~
%~~~~  Mean IR  ~~~~
hf1 = figure;
set(hf1,'Position',smallsq,'NextPlot','add')
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
set(gca,'xscale',axscale,'yscale',axscale)
xlim([axmin axmax])
ylim([axmin axmax])
xlabel(['Predicted IR ' USE_MEASURE])
ylabel(['Observed IR ' USE_MEASURE])


%~~~~~~~~~~~~~~~~~~~
%~~~~  Each IR  ~~~~
hf2 = figure;
set(hf2,'Position',smallsq,'NextPlot','add')
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
set(gca,'xscale',axscale,'yscale',axscale)
xlim([axmin axmax])
ylim([axmin axmax])
xlabel(['Predicted IR ' USE_MEASURE])
ylabel(['Observed IR ' USE_MEASURE])


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~  Diff vs baseline rate  ~~~~
hf3 = figure;
set(hf3,'Position',smallsq,'NextPlot','add')
subplot(1,5,1:4)
hold on
plot([0.01 10^2],[0 0],'Color',[0.7 0.7 0.7])
set(gca,'xscale','log')
xlim([axmin 10^2])
ylim([-5 5])
xlabel('Baseline FR')
ylabel(['Observed - Predicted IR ' USE_MEASURE])


% Set colors
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        
subjcolors = cmocean('phase',1+numel(unique({UnitData.Subject})));
Subjects   = unique({UnitData.Subject});

        
%% Select a datapoint (or a few?) to highlight

ex_subj = 'AAB_265058';
ex_sess = 'Jan17-AM';
ex_ch   = 63;
ex_clu  = 1230;

ex_Un   = find(strcmp({UnitData.Subject},ex_subj) & strcmp({UnitData.Session},ex_sess) & [UnitData.Channel]==ex_ch & [UnitData.Clu]==ex_clu);


%% Preallocate

N=0;
NIR=0;

Pdata = table;
Pdata.Subject = ' ';
Pdata.Session = ' ';
Pdata.Channel = nan;
Pdata.Clu     = nan;
Pdata.iUn     = nan;
Pdata.BaseFR  = nan;
Pdata.predIR  = nan;
Pdata.obsIR   = nan;
Pdata.stid    = nan;
Pdata.pval    = nan;




for iUn = 1:numel(UnitData)
    
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    shank       = UnitData(iUn).Shank;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    subjcol     = subjcolors(strcmp(subject,Subjects),:);
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
        
    
    % Load data files
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes and shift based on calculated integration time 
    
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 + spkshift));  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 + spkshift)');
        
    end
    
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    %%
    
    [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
        length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,20,'silence');
    
    % Find all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    
    if ~isfield('Info','artifact')
        [all_TDidx,Ntrials,~] = get_clean_trials(TrialData,[],dBSPL,LP);
    else
        [all_TDidx,Ntrials,~] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP);
    end
    
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
    FF      = nan(size(allStim'));
    FRtrials= cell(size(allStim));
    
    
    for istim = allStim'
        
        st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
        
        %%%  Skip ITI stimuli (?)
        ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
        TDidx = st_TDidx_ALL;%(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        t3 = TrialData.offset(TDidx);
        switch USE_MEASURE
            case 'FR'
                Duration = mode(diff([t2 t3],1,2));
            case 'FF'  % equalize analysis window across stimuli
                Duration = 1000;
        end
        t3 = t2 + Duration;
        if t3(end)>length(Stream_FRsmooth)
            t3 = t3(1:end-1);
            t2 = t2(1:end-1);
        end
        
        if exclOnset
            t2 = t2+150;
        end
        
        % Randomly permute the order of trials, for the simulation
        % below
        kt = randperm(length(t2));
        t2 = t2(kt);
        t3 = t3(kt);
        
        % Collect responses for each trial
        FR_resp = nan(numel(t2),1);
        
        for it = 1:numel(t2)
            
            switch units
                case 'Hz'
                    FR1 = mean(Stream_FRsmooth((t2(it)+1):t3(it)));
                    FR2 = sum(spiketimes>t2(it) & spiketimes<=t3(it)) / Duration * 1000;
                    %the outcome is identical either way
                    FR_resp(it,1) = FR1;
                case 'z'
                    FR_resp(it,1) = mean(Stream_zscore((t2(it)+1):t3(it)));
            end
            
        end %it
        
        
        % Calculate the mean FR for this stimulus
        FRmeans(istim==allStim')  = mean(FR_resp,1);
        FRvars(istim==allStim')   = var(FR_resp,1);
        
        % Calculate the mean FF for this stimulus
        FF(istim==allStim')       = bootstrap_for_FF(FR_resp,min(Ntrials(2:end)));
        % FF(istim)  = FRvars(istim) / FRmeans(istim);
        
        FRtrials{istim==allStim'} = FR_resp;
        
    end %istim
    
    
    %% ----  OBSERVED IR RESPONSE  ----
    
    switch USE_MEASURE
        case 'FR'
            IR_Observations  = FRmeans(7:end);
            IR_Obs_sems      = sqrt(FRvars(7:end)) ./ sqrt(Ntrials(7:end));
            
        case 'FF'
            IR_Observations  = FF(7:end);
            IR_Obs_sems      = zeros(size(IR_Observations));
    end
    
    
    
    %% ----  PREDICTED IR RESPONSE  ----
    
    % weighted average of periodic responses
    
    switch USE_MEASURE
        case 'FR'
            IR_Prediction   = sum( FRmeans(2:6) .* ((1./AMrates)/sum(1./AMrates)) );
            
            % Errorbar option #1: weighted sum of variances, then converted to overall SEM [--> var(a+b)=var(a)+var(b) ]
            IR_Pred_sem2   = sqrt( sum( FRvars(2:6) .* ((1./AMrates)/sum(1./AMrates)) ) ) / sqrt(length(AMrates));
            % Errorbar option #1: weighted sum of SEMs
            IR_Pred_sem    = sum( ( sqrt(FRvars(2:6))./sqrt(Ntrials(2:6)) ) .* ((1./AMrates)/sum(1./AMrates)) );
            
        case 'FF'
            IR_Prediction  = sum( FF(2:6) .* ((1./AMrates)/sum(1./AMrates)) );
            IR_Pred_sem    = 0;
    end
    
    
    
    %% Compare FR distributions
    
    if strcmp(USE_MEASURE,'FR')
        
        % Simulate IR trials for within unit statistics
        ntrs = cellfun(@length,FRtrials);
        nt = min(ntrs(ntrs>0));
        if nt<minTrs, keyboard, end
        
        FR_sim = zeros(nt,1);
        for it = 1:nt
            for ir = 1:5
                FR_sim(it) = FR_sim(it) + FRtrials{ir+1}(it) * ((1./AMrates(ir))/sum(1./AMrates));
            end
        end
        
        % Stats
        thisNIR = 0;
        pvals   = [];
        for iir = 7:length(FRtrials)
            if isempty(FRtrials{iir}), continue, end
            thisNIR = thisNIR+1;
            [~,pvals(thisNIR)] = ttest2(FR_sim,FRtrials{iir}(1:nt));
        end
        
    else
        pvals = nan(1,sum(allStim>6));
    end
    
    
    %% Save data to table
    
    for is = 1:sum(allStim>6)
        Pdata_addrow = { subject session channel clu iUn  UnitData(iUn).BaseFR  IR_Prediction  IR_Observations(is)  allStim(is+6)  pvals(is) };
        Pdata = [Pdata; Pdata_addrow];
    end
    
    
    %%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plotMTF && iUn==ex_Un                         % INDIVIDUAL UNIT MTF
        
        xvals = [1:5 8+(1:numel(allStim(7:end)))];
        switch USE_MEASURE
            case 'FR'
                y_vals = FRmeans';
                y_errs = FRmeans' + (FRvars.^0.5)'.*[-ones([length(FRmeans) 1]) ones([length(FRmeans) 1])]; %for sem instead of std ./ (Ntrials.^0.5)'
            case 'FF'
                y_vals = FF';
                y_errs = zeros(length(FF),2);
        end
        
        hmtf = figure; hold on
        set(hmtf,'Position',smallsq)
        
        % Plot baseline FR
        plot([0 max(xvals)+1],[UnitData(iUn).BaseFR UnitData(iUn).BaseFR],'--k','LineWidth',0.5)
        
        % Plot periodic stimuli
        for ir = 1:5
            plot(xvals(ir),y_vals(ir+1), 'o','MarkerSize',15,...
                'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
            plot([xvals(ir) xvals(ir)], y_errs(ir+1,:), ...
                '-','Color',colors(ir+1,:),'LineWidth',2)
        end
        
        % Plot IR prediction
        fill([6 max(xvals)+1 max(xvals)+1 6],IR_Prediction + IR_Pred_sem.*[-1 -1 1 1],0.8.*[1 1 1],'EdgeColor','none')
        plot([6 max(xvals)+1],[IR_Prediction IR_Prediction],'k','LineWidth',1)
        
        % Plot IR observed
        for ir = 6:length(xvals)
            plot(xvals(ir),y_vals(ir+1),'o','MarkerSize',15,...
                'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
            plot([xvals(ir) xvals(ir)],y_errs(ir+1,:),...
                '-','Color',colors(ir+1,:),'LineWidth',2)
        end
        
        % Finish formatting
        set(gca,'XTick',min(xvals):max(xvals),...
            'XTickLabel',[Info.stim_ID_key(2:6)' ' ' 'Linear Prediction' ' ' Info.stim_ID_key(allStim(7:end))' ],...
            'TickLabelInterpreter','none')
        xtickangle(45)
        
        xlim([-1+min(xvals) max(xvals)+1])
        ylim([axmin 15])
        ylabel(sprintf('%s response\n(%s)',USE_MEASURE, units))
        
        % Save MTF figure
        savedir = fullfile(fn.processed,subject,session,'Rasters',[subject '_' session '_' num2str(shank) '_' num2str(clu)]);
        if ~exist(savedir,'dir')
            keyboard
        end
        savename = sprintf('%s_%s_clu%i_FR-MTF_%s_exclOnset',...
            subject,session,clu,units);
%         print_eps_kp(hmtf,fullfile(savedir,savename),1)
%         print_svg_kp(hmtf,fullfile(savedir,savename),1)
        
    end
        
        
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                                                      POPULATION PLOT
    
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    % Add datapoint to plots
    % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    N = N+1;
    
    % Set color
    switch colorswitch
        case 'iBMF_FR'
            if isfield(UnitData(iUn),'iBMF_FR') && ~isempty(UnitData(iUn).iBMF_FR)
                plotcol = colors(1+UnitData(iUn).iBMF_FR,:);
            else
                plotcol = 'k';
            end
        case 'iBMF_VS'
            if isfield(UnitData(iUn),'iBMF_VS') && ~isempty(UnitData(iUn).iBMF_VS)
                plotcol = colors(1+UnitData(iUn).iBMF_VS,:);
            else
                plotcol = 'k';
            end
        case 'subject'
            plotcol = subjcol;
    end
    
    %~~~~~~~~~~~~~~~~~~~
    %~~~~  Mean IR  ~~~~
    figure(hf1); hold on
    
    % horizontal errorbars (prediction error)
    plot(IR_Prediction + IR_Pred_sem*[-1 1],[mean(IR_Observations) mean(IR_Observations)],'Color',plotcol)
    % vertical errorbars (observation error)
    plot([IR_Prediction IR_Prediction], mean(IR_Observations) + mean(IR_Obs_sems)*[-1 1],'Color',plotcol)
    
    % mean point
    ip=scatter(IR_Prediction,mean(IR_Observations),100,'o','LineWidth',2,'MarkerEdgeColor',plotcol,'MarkerFaceColor','none');
    if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
        ip.MarkerFaceColor = [32 129 255]./255;
    end
    if strcmp(USE_MEASURE,'FR') && any(pvals)<alfa
        ip.MarkerFaceColor = 'k';
    end
    
    %~~~~~~~~~~~~~~~~~~~
    %~~~~  Each IR  ~~~~
    figure(hf2); hold on
    
    for is = 1:sum(allStim>6)
        % horizontal errorbars (prediction error)
        plot(IR_Prediction + IR_Pred_sem*[-1 1],[IR_Observations(is) IR_Observations(is)],'Color',plotcol)
        % vertical errorbars (observation error)
        plot([IR_Prediction IR_Prediction], IR_Observations(is) + IR_Obs_sems(is)*[-1 1],'Color',plotcol)
        % mean point
        ip=scatter(IR_Prediction,IR_Observations(is),100,'o','LineWidth',2,'MarkerEdgeColor',plotcol,'MarkerFaceColor','none');
        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
            ip.MarkerFaceColor = [32 129 255]./255;
        end
        if strcmp(USE_MEASURE,'FR') && pvals(is)<alfa
            ip.MarkerFaceColor = 'k';
        end
        NIR = NIR+1;
    end %is
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~  Diff vs baseline rate  ~~~~
    figure(hf3); hold on
    
    for is = 1:sum(allStim>6)
        
        % mean point
        ip=scatter(UnitData(iUn).BaseFR, IR_Observations(is)-IR_Prediction, 200,'o','LineWidth',2,'MarkerEdgeColor',plotcol,'MarkerFaceColor','none');
        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
            ip.MarkerFaceColor = [32 129 255]./255;
        end
        if strcmp(USE_MEASURE,'FR') && pvals(is)<alfa
            ip.MarkerFaceColor = 'k';
        end
        
    end %is
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
end %iUn


%% Stats corrections?

Pdata(1,:)=[];

% [iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni([Pdata.pval],alfa);
% Pdata.iIR_sig_bonferroni     = iIR_sig_bonferroni';
% Pdata.iIR_sig_bonferroniholm = iIR_sig_bonferroniholm';


%% Finish figures

% hf1
figure(hf1); hold on
title(['Pred vs. Obs IR responses, avg IR, N = ' num2str(N) ' units'])

% hf2
figure(hf2); hold on
title(['Pred vs. Obs IR responses, each IR seq, N=' num2str(size(Pdata,1)) ' comparisons (' num2str(sum([Pdata.pval]<alfa)) ' sig), ' num2str(N) ' units'])
 
% hf3
xhist = -5:0.25:5;
diffhist = hist([Pdata.obsIR]-[Pdata.predIR],xhist);
yh = 10*ceil(max(diffhist)/10);

figure(hf3); 
subplot(1,6,6) 
hold on
fill(diffhist,xhist,'k','EdgeColor','none')
plot([yh yh],mean([Pdata.obsIR]-[Pdata.predIR],'omitnan') + std([Pdata.obsIR]-[Pdata.predIR],'omitnan').*[-1 1],'b-','LineWidth',2)
plot(yh,mean([Pdata.obsIR]-[Pdata.predIR],'omitnan'),'bd','MarkerSize',14,'MarkerEdgeColor','none','MarkerFaceColor','b')
xlim([0 yh])
ylim([-5 5])
set(gca,'ytick',[],'xtick',[])
box off


%% Save figures

savedir = fullfile(fn.figs,'PredObs');
if exclOnset
    savedir = fullfile(savedir,'exclOnset');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('PredObs_avgIR_%s',USE_MEASURE);
print_eps_kp(hf1,fullfile(savedir,savename))
print_svg_kp(hf1,fullfile(savedir,savename))

savename = sprintf('PredObs_eachIR_%s',USE_MEASURE);
print_eps_kp(hf2,fullfile(savedir,savename))
print_svg_kp(hf2,fullfile(savedir,savename))

savename = sprintf('PredObs_DiffsBaseFR_%s',USE_MEASURE);
print_eps_kp(hf3,fullfile(savedir,savename))
print_svg_kp(hf3,fullfile(savedir,savename))


save(fullfile(savedir,['Pdata_' USE_MEASURE]),'Pdata','-v7.3')



end


