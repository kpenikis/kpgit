function TRF_PredictedObserved_BS
%
%  TRF_PredictedObserved
%    version with splitting trials and bootstrapping 
%
%  KP, 2019-07
%


% close all

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
smth_win    = 10;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Irr         = [6 7];
Pdc         = [1 2 3 4 5]; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TLAGS       = round(logspace(log10(20),log10(500),10));
% TLAGS       = [];
d_tlag      = 100;
TLAGS       = sort([TLAGS d_tlag]);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
exclOnset   = 0; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%% Figure settings

savedir = fullfile(fn.figs,'TRF');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
subjcolors = cmocean('phase',1+numel(unique({UnitData.Subject})));
Subjects   = unique({UnitData.Subject});

AMrates   = [2 4 8 16 32];
Durations = [1500 1000 1000 1000 1000 1000 1937 1937 1000];

rng('shuffle')

% Preallocate
Results = table;
Results.Subject  = ' ';
Results.Session  = ' ';
Results.Channel  = nan;
Results.Clu      = nan;
Results.iUn      = nan;
Results.C100_Pdc = nan;
Results.C100_Irr = nan;
Results.p_val    = nan;
Results.corrBest = nan;
Results.tlagBest = nan;
Results.BaseFR   = nan;
Results.drivFR   = nan;
Results.maxVS    = nan;
Results.maxVSrt  = nan;
Results.BMFrt    = nan;

C100_Pdc =[];
C100_Irr =[];
C500_Pdc =[];
C500_Irr =[];
Cmax_Pdc =[];
Cmax_Irr =[];
Cmax_ALL =[];
Imax_ALL =[];
STRF_corr=[];

for iUn = 1:numel(UnitData)
    
%     if isempty([UnitData(iUn).iBMF_FR  UnitData(iUn).iBMF_VS])
%         clear Phase0 Clusters Spikes
%         continue
%     end
    
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
    
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 ||  ~exist('Phase0','var') 
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info 
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1 ||  ( ~exist('Spikes','var') || ~exist('Clusters','var') )
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    
    % Get spiketimes and shift based on calculated integration time 
    
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift));  %ms
%         spiketimes = unique(round(  (Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift) / 1000 * FSfs )); %samples
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
        
    end
    
    % If shuffling spiketimes
%     spiketimes = sort(randi( length(SpoutStream), size(spiketimes) ));
    
    % Convert to samples
%     TrialData.onset  = round([TrialData.onset]/1000*FSfs);
%     TrialData.offset = round([TrialData.offset]/1000*FSfs);
%     TLAGS = round(TLAGS/1000*FSfs);
    
%     fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    
    %%
    
    [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
        length(SoundStream),TrialData.onset(1),TrialData.offset(1),[],smth_win,'silence');
    
    % Find all stimuli presented with these parameters, given a
    % sufficient number of trials without diruptive artifact
    % while the animal was drinking
    
    [all_TDidx,Ntrials,~,allStim] = get_clean_trials(TrialData,Info.artifact(channel).trials,dBSPL,LP,0);
    
    Ntrials(1) = 200; %set high so Warn stimulus doesn't interfere
    
    % N trials of Stim and PSTH per group
    nTrGrp = min(Ntrials(2:6)); %floor(min(Ntrials(2:end))/2);
    if nTrGrp<16 && iUn~=87
        keyboard
    end
    
    if any(Ntrials(7:end) < nTrGrp*2)
        thisIR = allStim(Ntrials(7:end) < nTrGrp*2) + 6;
        all_TDidx(TrialData.trID(all_TDidx)==thisIR)  = [];
        Ntrials(allStim==thisIR) = [];
        allStim(allStim==thisIR) = [];
    end
    
    
    
    %% 
    
    STIM_cell = cell(7,1);
    PSTH_cell = cell(7,1);
    
    for stid = allStim(allStim>1)'
        
        ist = find(stid==allStim);
        
        st_TDidx_ALL = all_TDidx(TrialData.trID(all_TDidx)==stid);
        
        %%%  Skip ITI stimuli (?)
        ITIflag = 0;%unique(TrialData.ITIflag(st_TDidx_ALL));
        TDidx = st_TDidx_ALL(TrialData.ITIflag(st_TDidx_ALL) == ITIflag(1));
        
        % Get timestamps of onsets and offsets
        clear t2 t3 Duration t_win
        t2 = TrialData.onset(TDidx);
        Duration = Durations(stid); %mode(diff([t2 t3],1,2));
        t3 = t2 + Duration;
        if t3(end)>length(Stream_FRsmooth)
            t3 = t3(1:end-1);
            t2 = t2(1:end-1);
        end
        
        if exclOnset
            t2 = t2+150;
        end
        
        % Permute order of trials
%         kt = randperm(length(t2));
%         t2 = t2(kt);
%         t3 = t3(kt);
        
        
        % Collect responses for each trial
        stim_i   = nan(  numel(t2), Duration+1);
        psth_rw  = zeros(numel(t2), Duration+1);
        psth_sm  = nan(  numel(t2), Duration+1);
        
        for it = 1:numel(t2)
            
            stim_i(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
                        
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) +1 ;
            
            psth_rw(it,sp) = 1000;  % raw, compare to smoothed 
            psth_sm(it,:) = Stream_FRsmooth(1, t2(it) : t3(it) );
            
        end %it
        
        STIM_cell{ist-1} = exp( stim_i );
        PSTH_cell{ist-1} = psth_sm;
        
    end %istim
    
    
    % Normalize envelope across stimuli
    stimMin = min(cell2mat(cellfun(@(x) min(min(x)),STIM_cell,'UniformOutput',false)));
    foo = cellfun(@(x) x-stimMin, STIM_cell,'UniformOutput',false);
    %     STIM = STIM/max(max(STIM));
    STIM_cell = foo; clear foo
    
    
    % OPTIONAL: Subtract spontaneous FR from PSTH
%     foo = cellfun(@(x) max(x-UnitData(iUn).BaseFR,0), PSTH_cell, 'UniformOutput',false);
%     PSTH_cell = foo; clear foo
    
    

    %%  ------  STRF analysis  ------ 
    
    nBS = 100;
    
    STRF_Pdc      = nan(numel(TLAGS),max(TLAGS),nBS);
    STRF_Irr      = nan(numel(TLAGS),max(TLAGS),nBS);
    coinc_Pdc_Pdc = nan(numel(TLAGS),nBS);
    coinc_Pdc_Irr = nan(numel(TLAGS),nBS);
    coinc_Irr_Pdc = nan(numel(TLAGS),nBS);
    coinc_Irr_Irr = nan(numel(TLAGS),nBS);
    
    for iBS = 1:nBS
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Concatenate PERIODIC data
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        psth_Pdc1 = []; stim_Pdc1 = [];
        
        for st = Pdc
            
            allTrs = randperm(size(STIM_cell{st},1));
            if isempty(allTrs), keyboard, end
            
            % Concatenate trials
%             psth_Pdc1  = [psth_Pdc1 reshape( PSTH_cell{st}(allTrs(1:nTrGrp),:)' , 1, numel(PSTH_cell{st}(allTrs(1:nTrGrp),:)) ) ];
%             stim_Pdc1  = [stim_Pdc1 reshape( STIM_cell{st}(allTrs(1:nTrGrp),:)' , 1, numel(STIM_cell{st}(allTrs(1:nTrGrp),:)) ) ];
            
            % Average trials
            psth_Pdc1  = [psth_Pdc1 mean(PSTH_cell{st}(allTrs(1:nTrGrp),:),1) ];
            stim_Pdc1  = [stim_Pdc1 mean(STIM_cell{st}(allTrs(1:nTrGrp),:),1) ];
            
        end
        if ~any(find(psth_Pdc1)) %|| ~any(find(psth_Pdc2))
            continue
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Concatenate IRREGULAR data
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        psth_Irr1 = []; stim_Irr1 = [];
        psth_Irr2 = []; stim_Irr2 = [];
        
        for st = Irr
            
            allTrs = randperm(size(STIM_cell{st},1));
            if isempty(allTrs), continue, end
            
            % Concatenate trials
%             psth_Irr1  = [psth_Irr1 reshape( PSTH_cell{st}(allTrs(1:nTrGrp),:)', 1, numel(PSTH_cell{st}(allTrs(1:nTrGrp),:)) ) ];
%             stim_Irr1  = [stim_Irr1 reshape( STIM_cell{st}(allTrs(1:nTrGrp),:)', 1, numel(STIM_cell{st}(allTrs(1:nTrGrp),:)) ) ];
%             
%             psth_Irr2  = [psth_Irr2 reshape( PSTH_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:)', 1, numel(PSTH_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:)) ) ];
%             stim_Irr2  = [stim_Irr2 reshape( STIM_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:)', 1, numel(STIM_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:)) ) ];
            
            % Average trials
            psth_Irr1  = [psth_Irr1 mean(PSTH_cell{st}(allTrs(1:nTrGrp),:),1) ];
            stim_Irr1  = [stim_Irr1 mean(STIM_cell{st}(allTrs(1:nTrGrp),:),1) ];
                        
            psth_Irr2  = [psth_Irr2 mean(PSTH_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
            stim_Irr2  = [stim_Irr2 mean(STIM_cell{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
        
        end
        if ~any(find(psth_Irr1)) || ~any(find(psth_Irr2))
            continue
        end
        
        
        % For each tlag, calculate STRF, predict response, and compare to observed
        
        for it = 1:numel(TLAGS)
            
            tlag  = TLAGS(it);
            
            %......................................
            % Estimate STRF from PERIODIC stimuli
            %......................................
            pred_response=[]; rP_norm=[]; rO_norm=[];
            STRF_Pdc(it,1:tlag,iBS) = sub_revCOR_AM(psth_Pdc1,stim_Pdc1,tlag);
            
            %......................................
            % Estimate STRF from IRREGULAR stimuli
            %......................................
            pred_response=[]; rP_norm=[]; rO_norm=[];
            STRF_Irr(it,1:tlag,iBS) = sub_revCOR_AM(psth_Irr1,stim_Irr1,tlag);
            
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Predict response to IRREGULAR stimuli
            %......................................
            % From Pdc STRF
            %.......................................
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_Irr2,STRF_Pdc(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_Irr2/norm(psth_Irr2);
            
            % Correlation between predicted and observed responses
            coinc_Pdc_Irr(it,iBS) = max( xcorr( rP_norm, rO_norm, 0) );
            
            
            %.......................................
            % From Pdc STRF
            %.......................................
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_Irr2,STRF_Irr(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_Irr2/norm(psth_Irr2);
            
            % Correlation between predicted and observed responses
            coinc_Irr_Irr(it,iBS) = xcorr( rP_norm, rO_norm, 0);
            
            
        end %tlags
    end %iBS
    
    if mode(sum(isnan(coinc_Irr_Irr),2)/size(coinc_Irr_Irr,2)) > 0.5
        keyboard
        coinc_Irr_Irr = nan(size(coinc_Irr_Irr));
    end
    
    % Best tlag from Irr Irr match
    [c_max,i_max] = max(mean(coinc_Irr_Irr,2,'omitnan'));
    
    C100_Pdc = [C100_Pdc; mean(coinc_Pdc_Irr(TLAGS==100,:),'omitnan')];
    C100_Irr = [C100_Irr; mean(coinc_Irr_Irr(TLAGS==100,:),'omitnan')];
    C500_Pdc = [C500_Pdc; mean(coinc_Pdc_Irr(TLAGS==500,:),'omitnan')];
    C500_Irr = [C500_Irr; mean(coinc_Irr_Irr(TLAGS==500,:),'omitnan')];
    Cmax_Pdc = [Cmax_Pdc; mean(coinc_Pdc_Irr(i_max,:),'omitnan')];
    Cmax_Irr = [Cmax_Irr; mean(coinc_Irr_Irr(i_max,:),'omitnan')];
    
    Cmax_ALL = [Cmax_ALL; c_max];
    Imax_ALL = [Imax_ALL; i_max];
    
    p_tt = [];
    [~,p_tt] = ttest2(coinc_Pdc_Irr(TLAGS==100,:),coinc_Irr_Irr(TLAGS==100,:));
    
    
    %% Save data to table
    
    VSmax = max( UnitData(iUn).VSdata_spk( 1, UnitData(iUn).VSdata_spk(3,:)<0.05 ) );
    if isempty(VSmax)
        VSmax=nan;
        VSrt=nan;
    else
        VSrt = AMrates(find(UnitData(iUn).VSdata_spk(1,:)==VSmax) - 1);
    end
    BMFrt = AMrates(UnitData(iUn).iBMF_FR);
    if isempty(BMFrt)
        BMFrt=nan;
    end
    
    % Override any data for this unit
    if iUn==87
        for ii=6:10
            Results_addrow{ii}=nan;
        end
    end
    
    % Add to Results table
    Results_addrow = { subject session channel clu  iUn   mean(coinc_Pdc_Irr(TLAGS==100,:),'omitnan')   mean(coinc_Irr_Irr(TLAGS==100,:),'omitnan')  p_tt  c_max  i_max   UnitData(iUn).BaseFR  mean(psth_Pdc1)  VSmax  VSrt  BMFrt  };
    Results = [Results; Results_addrow];
    
end %iUn

Results(1,:)=[];



%% Plot results
% train Pdc predict Irr   vs.   train Irr predict Irr

% Load results from Pred vs Obs avg FR
loaddir = fullfile(fn.figs,'PredObs');
Pdata_FR = load(fullfile(loaddir,'Pdata_FR'));
Pdata_FR = Pdata_FR.Pdata;

% Set color according to difference from FR PredObs analysis
unNL = nan(numel(UnitData),1);
for iu = 1:numel(UnitData)
    unNL(iu) = mean(Pdata_FR(Pdata_FR.iUn==iu,:).obsIR - Pdata_FR(Pdata_FR.iUn==iu,:).predIR);
end

axmax = 1;


%.......................................
% Focus on fixed tlag of 100
%.......................................

hf100=figure;
set(hf100,'Position',fullscreen,'NextPlot','add')

% Scatter plot
hs(1)=subplot(2,3,1);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on
try
    scatter(C100_Pdc,C100_Irr,100,abs(unNL),'LineWidth',1)
    cmocean('thermal')
catch
    scatter(C100_Pdc,C100_Irr,100,'k','LineWidth',1)
end
set(gca,'Color','none')
title('tlag=100')
xlabel('STRF from Pdc')
ylabel('STRF from Irr')

% Histogram of x axis
hs(4)=subplot(2,3,4);
histogram(C100_Pdc,0:0.05:1,'FaceColor','k','EdgeColor','none','FaceAlpha',1,'Normalization','probability')
xlabel('Corr from Pdc STRF')
title('Histogram of x axis, tlag=100')
axis square; box off
set(gca,'Color','none')
ylim([0 0.5])

% Histogram of y axis
hs(5)=subplot(2,3,5);
histogram(C100_Irr,0:0.05:1,'FaceColor','k','EdgeColor','none','FaceAlpha',1,'Normalization','probability')
xlabel('Corr from Irr STRF')
title('Histogram of y axis, tlag=100')
axis square; box off
set(gca,'Color','none')
ylim([0 0.5])

% Histogram of difference
hs(6)=subplot(2,3,6);
histogram(C100_Pdc-C100_Irr,-0.5:0.05:0.5,'FaceColor','k','EdgeColor','none','FaceAlpha',1,'Normalization','probability')
xlabel('Corr Irr - Corr Pdc')
title('Histogram of diagonal, tlag=100')
axis square; box off
set(gca,'Color','none')


% Use extra subplots to display other tlags

hs(2)=subplot(2,3,2);
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
hold on
scatter(C500_Pdc,C500_Irr,100,'b','LineWidth',1)
title('tlag=500')
xlabel('STRF from Pdc')
ylabel('STRF from Irr')
axis square; box on
set(gca,'Color','none')

hs(3)=subplot(2,3,3);  %shows that the bias isn't that bad
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
hold on
scatter(Cmax_Pdc,Cmax_Irr,100,'r','LineWidth',1)
title('Best tlag for Irr')
xlabel('STRF from Pdc')
ylabel('STRF from Irr')
axis square; box on
set(gca,'Color','none')


suptitle('Coincidence values for predicting IR response, from STRFs from Pdc vs. Irr stimulus')


savename = sprintf('PredictIR');
print_eps_kp(hf100,fullfile(savedir,savename));


% Also save Pdata table
save(fullfile(savedir,'Pdata_TRF'),'Results','-v7.3')



keyboard


end


