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


% close all
global fn AMrates rateVec_AC rateVec_DB 

%!!!!!!!!!!!!!!!!!
exclOnset   = 1; keyboard
%!!!!!!!!!!!!!!!!!
plotMTF     = 1;
%!!!!!!!!!!!!!!!!!
minTrs      = 10;
%!!!!!!!!!!!!!!!!!
colorswitch = 'iBMF_FR'; 'iBMF_VS'; 'subject';
%!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


if nargin<1
    USE_MEASURE = 'FR'; 'TrV'; 'FF';
end

switch USE_MEASURE
    case 'FR'
        units = 'Hz';
        switch units
            case 'Hz'
                axmin = 0.01;
                axmax = 10*ceil(max([UnitData.BaseFR])/10);
            case 'z'
                axmin = -1;
                axmax = 1.5;
        end
    case {'TrV' 'FF'}
        units = ' ';
        axmin = 0.01;
        axmax = 4;
end



%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');
vsmallsq = [1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/4];
smallsq  = [1 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3];

% Prepare the figure
hf1 = figure;
set(hf1,'Position',smallsq,'NextPlot','add')
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
set(gca,'xscale','log','yscale','log')
xlim([axmin axmax])
ylim([axmin axmax])
xlabel('Predicted')
ylabel('Observed')

% Prepare another figure, for individual sequences' observed responses
hf2 = figure;
set(hf2,'Position',smallsq,'NextPlot','add')
hold on
plot([axmin axmax],[axmin axmax],'Color',[0.7 0.7 0.7])
axis square
set(gca,'xscale','log','yscale','log')
xlim([axmin axmax])
ylim([axmin axmax])
xlabel('Predicted IR response')
ylabel('Observed IR response')

% Prepare a third figure, Observed-Predicted as a fct of BaseFR
hf3 = figure;
set(hf3,'Position',smallsq,'NextPlot','add')
subplot(1,5,1:4)
hold on
plot([axmin axmax],[0 0],'Color',[0.7 0.7 0.7])
set(gca,'xscale','log')
xlim([axmin axmax])
ylim([-5 5])
xlabel('Baseline FR')
ylabel('Observed - Predicted IR response')

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
Pdata.IRstim  = nan;
Pdata.BaseFR  = nan;
Pdata.pval    = nan;
Pdata.predIR  = nan;
Pdata.obsIR   = nan;



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
    subjcol     = [1 1 1];
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
    
    spkshift    = UnitData(iUn).IntTime_spk;
    
    
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
    FRtrials= cell(size(allStim));
    
    for istim = allStim'
        
        st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==istim);
        
        %%%  Skip ITI stimuli (?)
        ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
        
        for is = 1:numel(ITIflag)
            TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(is));
            
            % Get timestamps of onsets and offsets
            clear t2 t3 Duration t_win
            t2 = TrialData.onset(TDidx);
            t3 = TrialData.offset(TDidx);
            Duration = mode(diff([t2 t3],1,2));
            t3 = t2 + Duration;
            
            if exclOnset
                t2 = t2+200;
            end
            
            % Randomly permute the order of trials, for the simulation
            % below
            kt = randperm(length(t2));
            t2 = t2(kt);
            t3 = t3(kt);
            
            % Now collect FR responses
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
            
            % -- OBSERVED --
            % Calculate the mean normalized FR for this stimulus
            FRmeans(istim) = mean(FR_resp,1);
            FRvars(istim)  = var(FR_resp,1);
            
            FRtrials{istim} = FR_resp;
            
        end %iti filter
    end %istim
    
    
    
    
    %% ----  PREDICTED IR RESPONSE  ----
    
    % Calculate prediction for IR response, based on
    % weighted average of periodic responses
    IR_Prediction = sum( FRmeans(2:6) .* ((1./AMrates)/sum(1./AMrates)) );
    
    % add error bars
    %weighted sum of variances, then converted to overall SEM [--> var(a+b)=var(a)+var(b) ]
    IR_Pred_sem2 = sqrt( sum( FRvars(2:6) .* ((1./AMrates)/sum(1./AMrates)) ) ) / sqrt(length(AMrates));
    %weighted sum of SEMs
    IR_Pred_sem  = sum( ( sqrt(FRvars(2:6))./sqrt(Ntrials(2:6)) ) .* ((1./AMrates)/sum(1./AMrates)) );
    
    
    
    %% Simulate IR trials for within unit statistics
    
    ntrs = cellfun(@length,FRtrials);
    nt = min(ntrs(ntrs>0));
    if nt<minTrs, keyboard, end
    
    FR_sim = zeros(nt,1);
    for it = 1:nt
        for ir = 1:5
            FR_sim(it) = FR_sim(it) + FRtrials{ir+1}(it) * ((1./AMrates(ir))/sum(1./AMrates));
        end
    end
    
    if ~isempty(FRtrials{7})
        NIR = NIR+1;
        
        [~,pvals(NIR)] = ttest2(FR_sim,FRtrials{7}(1:nt));
        
        % Save data to table
        Pdata_addrow = { subject session channel clu 7 UnitData(iUn).BaseFR pvals(NIR) IR_Prediction mean(FRtrials{7})  };
        Pdata = [Pdata; Pdata_addrow];
    end
    if numel(FRtrials)>7
        NIR = NIR+1;
        
        [~,pvals(NIR)] = ttest2(FR_sim,FRtrials{8}(1:nt));
        
        % Save data to table
        Pdata_addrow = { subject session channel clu 8 UnitData(iUn).BaseFR pvals(NIR) IR_Prediction mean(FRtrials{8})  };
        Pdata = [Pdata; Pdata_addrow];
    end
    
    
    
    
    %%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if plotMTF && iUn==ex_Un                     % INDIVIDUAL UNIT ZFR MTF
        
        xvec = [1:5 8+(1:numel(allStim(7:end)))];
        msize = Ntrials./5;
        
        hmtf = figure; hold on
        set(hmtf,'Position',smallsq)
        
        % Plot baseline FR
        plot([0 max(xvec)+1],[UnitData(iUn).BaseFR UnitData(iUn).BaseFR],'--k','LineWidth',0.5)
        
        % Plot periodic stimuli
        for ir = 1:5
            plot(xvec(ir),FRmeans(ir+1),'o','MarkerSize',15,...
                'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
            plot([xvec(ir) xvec(ir)],[FRmeans(ir+1)-(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1))) FRmeans(ir+1)+(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1)))],...
                '-','Color',colors(ir+1,:),'LineWidth',2)
        end
        
        % Plot IR prediction
        fill([6 max(xvec)+1 max(xvec)+1 6],IR_Prediction + IR_Pred_sem.*[-1 -1 1 1],0.8.*[1 1 1],'EdgeColor','none')
        plot([6 max(xvec)+1],[IR_Prediction IR_Prediction],'k','LineWidth',1)
%         plot(7,IR_Prediction,'o','MarkerSize',15,'LineWidth',2,'Color',[32 129 255]./255)
%         plot([7 7],[IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],'-','LineWidth',2,'Color',[32 129 255]./255)
        
        % Plot IR observed
        for ir = 6:length(xvec)
            plot(xvec(ir),FRmeans(ir+1),'o','MarkerSize',15,...
                'MarkerFaceColor',colors(ir+1,:),'MarkerEdgeColor','none')
            plot([xvec(ir) xvec(ir)],[FRmeans(ir+1)-(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1))) FRmeans(ir+1)+(sqrt(FRvars(ir+1))/sqrt(Ntrials(ir+1)))],...
                '-','Color',colors(ir+1,:),'LineWidth',2)
        end
        
        % Finish formatting
        set(gca,'XTick',min(xvec):max(xvec),...
            'XTickLabel',[Info.stim_ID_key(2:6)' ' ' 'Linear Prediction' ' ' Info.stim_ID_key(allStim(7:end))' ],...
            'TickLabelInterpreter','none')
        xtickangle(45)
        
        xlim([-1+min(xvec) max(xvec)+1])
        ylim([axmin 15])
        ylabel(sprintf('FR response\n(%s)', units))
        
        % Save MTF figure
        savedir = fullfile(fn.processed,subject,session,'Rasters',[subject '_' session '_' num2str(shank) '_' num2str(clu)]);
        if ~exist(savedir,'dir')
            keyboard
        end
        savename = sprintf('%s_%s_clu%i_FR-MTF_%s_exclOnset',...
            subject,session,clu,units);
        print_eps_kp(hmtf,fullfile(savedir,savename),1)
%         print_svg_kp(hmtf,fullfile(savedir,savename),1)
        
    end
        
        
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %                                            POPULATION ZFR COMPARISON
    
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
    
    
    figure(hf1); hold on
    
    % horizontal errorbars (prediction error)
    plot([IR_Prediction-IR_Pred_sem IR_Prediction+IR_Pred_sem],[mean(FRmeans(7:end)) mean(FRmeans(7:end))],'Color',plotcol)
    % vertical errorbars (observation error)
    plot([IR_Prediction IR_Prediction],[mean(FRmeans(7:end))-mean(sqrt(FRvars(7:end))./sqrt(Ntrials(7:end)))...
        mean(FRmeans(7:end))+mean(sqrt(FRvars(7:end))./sqrt(Ntrials(7:end)))],'Color',plotcol)
    %                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10))) mean(zFRmeans(7:10))+sqrt(mean(zFRvars(7:10)))./sqrt(sum(blocks_N(7:10)))],plotcol)
    %                     plot([IR_Prediction IR_Prediction],[mean(zFRmeans(7:10))-sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4) mean(zFRmeans(7:10))+sqrt((1/4)*sum(zFRvars(7:10)))./sqrt(4)],plotcol)
    % mean point
    ip=scatter(IR_Prediction,mean(FRmeans(7:end)),150,'o','MarkerFaceAlpha',0.45,'MarkerFaceColor',plotcol,'MarkerEdgeColor','none');
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
        ip=scatter(IR_Prediction,FRmeans(is),100,'o','MarkerFaceAlpha',0.45,'MarkerFaceColor',plotcol,'MarkerEdgeColor','none');
        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
            ip.MarkerFaceColor = [32 129 255]./255;
            ip.MarkerFaceAlpha = 0.65;
        end
        
    end %is
    
    
    figure(hf3); hold on
    
    for is = 7:length(allStim)
        
        % vertical errorbars (observation error)
        %             plot([UnitData(iUn).BaseFR UnitData(iUn).BaseFR], FRmeans(is)-IR_Prediction + (sqrt(FRvars(is))./sqrt(Ntrials(is))).*[-1 1],'Color',plotcol)
        % mean point
        ip=scatter(UnitData(iUn).BaseFR, FRmeans(is)-IR_Prediction, 200,'o','MarkerFaceAlpha',0.45,'MarkerFaceColor',plotcol,'MarkerEdgeColor','none');
        if strcmp(subject,ex_subj) && strcmp(session,ex_sess) && channel==ex_ch && clu==ex_clu
            ip.MarkerFaceColor = [32 129 255]./255;
            ip.MarkerFaceAlpha = 0.65;
        end
        
    end %is
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
end %iUn


%% Stats
[iIR_sig_bonferroni,iIR_sig_bonferroniholm] = checkSignificance_bonferroni(pvals,0.05);

Pdata(1,:)=[];

Pdata.iIR_sig_bonferroni     = iIR_sig_bonferroni';
Pdata.iIR_sig_bonferroniholm = iIR_sig_bonferroniholm';


%% Finish figures

% hf1
figure(hf1); hold on
title(['Pred vs. Obs IR responses, avg IR, N = ' num2str(N) ' units'])

% hf2
figure(hf2); hold on
plot(Pdata.predIR(find(iIR_sig_bonferroni))',Pdata.obsIR(find(iIR_sig_bonferroni))','sg')
title(['Pred vs. Obs IR responses, each IR seq, N=' num2str(size(Pdata,1)) ', ' num2str(N) ' units'])

% hf3
xhist = -5:0.25:5;
diffhist = hist([Pdata.obsIR]-[Pdata.predIR],xhist);
yh = 10*ceil(max(diffhist)/10);

figure(hf3); 
subplot(1,6,6) 
hold on
fill(diffhist,xhist,'k','EdgeColor','none')
plot([yh yh],mean([Pdata.obsIR]-[Pdata.predIR]) + std([Pdata.obsIR]-[Pdata.predIR]).*[-1 1],'b-','LineWidth',3)
plot(yh,mean([Pdata.obsIR]-[Pdata.predIR]),'bd','MarkerSize',14,'MarkerEdgeColor','none','MarkerFaceColor','b')
xlim([0 yh])
ylim([-5 5])
set(gca,'ytick',[],'xtick',[])
box off


%% Save figures

savedir = fullfile(fn.figs,'PredObs','exclOnset');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = 'PredObs_avgIR';
print_eps_kp(hf1,fullfile(savedir,savename))
print_svg_kp(hf1,fullfile(savedir,savename))

savename = 'PredObs_eachIR';
print_eps_kp(hf2,fullfile(savedir,savename))
print_svg_kp(hf2,fullfile(savedir,savename))

savename = 'PredObs_DiffsBaseFR';
print_eps_kp(hf3,fullfile(savedir,savename))
print_svg_kp(hf3,fullfile(savedir,savename))


save(fullfile(savedir,'Pdata'),'Pdata','-v7.3')



end


