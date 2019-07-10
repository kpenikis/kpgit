function TRF_PredictedObserved(FullSound,FSfs)
%
%  TRF_PredictedObserved
%
%  KP, 2019-06
%

% for significance calculation must measure each IR stimulus separately
% (add in fig 2 and run again)


% close all

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
smth_win    = 10;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Irr          = [6 7];
Pdc          = [1 2 3 4 5]; 
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

AMrates   = [2 4 8 16 32];
Durations = [1500 1000 1000 1000 1000 1000 1937 1937 1000];

scrsz = get(0,'ScreenSize');  %[left bottom width height]
widescreen   = [1 scrsz(4) scrsz(3) scrsz(4)/3*2];
set(0,'DefaultAxesFontSize',12)

% Preallocate
Pdata = table;
Pdata.Subject  = ' ';
Pdata.Session  = ' ';
Pdata.Channel  = nan;
Pdata.Clu      = nan;
Pdata.iUn      = nan;
Pdata.Corrs    = cell(1,1);
Pdata.corr0    = nan;
Pdata.corrBest = nan;
Pdata.tlagBest = nan;
Pdata.Range    = nan;
Pdata.BaseFR   = nan;
Pdata.drivFR   = nan;
Pdata.maxVS    = nan;
Pdata.maxVSrt  = nan;
Pdata.BMFrt    = nan;

rng('shuffle')

Cmax_ALL=[];
Imax_ALL=[];
STRF_corr=[];

for iUn = 124:numel(UnitData)
    
%     close all
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    shank       = UnitData(iUn).Shank;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
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
    nTrGrp = floor(min(Ntrials(2:end))/2);
    
    if any(Ntrials(2:6) < nTrGrp)
        continue
    end
    
    if sum(Ntrials < nTrGrp)==1
        all_TDidx(TrialData.trID(all_TDidx)==allStim(Ntrials<nTrGrp))  = [];
        allStim(Ntrials<nTrGrp)  = [];
        Ntrials(Ntrials<nTrGrp) = [];
    elseif  sum(Ntrials < nTrGrp)>1
        keyboard
    end
    
    
    %% 
    
    STIM1 = nan(7,2000);
    PSTH1 = nan(7,2000);
    STIM2 = nan(7,1e7);
    PSTH2 = nan(7,1e7);
    
    STIM_sp = nan(14,2000);
    PSTH_sp = nan(14,2000);
    
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
%         stim_i   = zeros(numel(t2), Duration+1);
%         psth_rw  = zeros(numel(t2), Duration+1);
%         psth_sm  = nan(  numel(t2), Duration+1);
        
        stimcat  = [];
        psthcat  = [];
        
        stim_i   = zeros(nTrGrp, Duration+1);
        psth_sm  = nan(  nTrGrp, Duration+1);
        
        grp=0;
        for it = 1:numel(t2)
            
            stim_i(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
%             stimcat = [stimcat SoundStream(1, t2(it) : t3(it) )...
%                 ./ max(SoundStream(1, t2(it) : t3(it) )) ];
%             stimcat = [stimcat FullSound(1, t2(it) : t3(it) ) ];
            
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) +1 ;
            
            psth_rw(it,sp) = 1000;  % raw, compare to smoothed 
            psth_sm(it,:) = Stream_FRsmooth(1, t2(it) : t3(it) );
            
            psthcat = [psthcat Stream_FRsmooth(1, t2(it) : t3(it) )];
            
            if (it==nTrGrp) || (it==nTrGrp*2)
                STIM_sp( (ist-1)+6*grp, 1:Duration+1) = mean(stim_i,1,'omitnan');
                PSTH_sp( (ist-1)+6*grp, 1:Duration+1) = mean(psth_sm,1,'omitnan');
                grp = grp+1;
            end
            
        end %it
        
        STIM1(ist-1,1:Duration+1) = mean(stim_i,1,'omitnan');
        PSTH1(ist-1,1:Duration+1) = mean(psth_sm,1,'omitnan');
        
        STIM2(ist-1,1:length(stimcat)) = stimcat;
        PSTH2(ist-1,1:length(psthcat)) = psthcat;
        
    end %istim
    
%     STIM = STIM1;
%     PSTH = PSTH1;
    
    % Normalize envelope across all stim
    STIM = exp( STIM1 );
    STIM = (STIM-min(min(STIM)));
%     STIM = STIM/max(max(STIM));
    
    % Baseline-subtracted, rectified FR
    PSTH = PSTH1;
    baseFR = UnitData(iUn).BaseFR;
%     drivFR = mean(1000*sum(spiketimes>TrialData.onset(all_TDidx) & spiketimes<TrialData.offset(all_TDidx),2)./(TrialData.offset(all_TDidx)-TrialData.onset(all_TDidx)));
%     PSTH   = max( PSTH-baseFR , 0);
    
    




    %%  ------  STRF analysis  ------ 
    
    % Concatenate training data
    psth_Pdc = []; stim_Pdc = [];
%     psth_Pdc = zeros(1,500);
%     stim_Pdc = zeros(1,500);
    for t = Pdc  
        s = sum(isnan(STIM(t,:)));
        psth_Pdc  = [psth_Pdc PSTH(t,1:(end-s)) ];
        stim_Pdc  = [stim_Pdc STIM(t,1:(end-s)) ];
    end
    
%     FFTstim_Pdc = fft(stim_Pdc);
%     SpecStim = spectrogram(stim_Pdc,round(FSfs/200),round(FSfs/500),logspace(log10(1),log10(round(FSfs/2)),10),FSfs,'yaxis');
    
    % Concatenate testing data
    psth_Irr = []; stim_Irr = [];
%     psth_Irr = zeros(1,500);
%     stim_Irr = zeros(1,500);
    for t = Irr
        s = sum(isnan(STIM(t,:)));
        psth_Irr  = [psth_Irr PSTH(t,1:(end-s)) ];
        stim_Irr  = [stim_Irr STIM(t,1:(end-s)) ];
    end
    
    
    % For each tlag, calculate STRF, predict response, and compare to observed
    
    STRF_Pdc      = nan(numel(TLAGS),max(TLAGS));
    STRF_Irr      = nan(numel(TLAGS),max(TLAGS));
    coinc_Pdc_Pdc = nan(numel(TLAGS),1);
    coinc_Pdc_Irr = nan(numel(TLAGS),1);
    coinc_Irr_Pdc = nan(numel(TLAGS),1);
    coinc_Irr_Irr = nan(numel(TLAGS),1);
%     STRF_test     = nan(numel(TLAGS),max(TLAGS));
%     coinc_test    = nan(numel(TLAGS),1);
    
    for it = 1:numel(TLAGS)
        
        tlag  = TLAGS(it);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Estimate STRF from PERIODIC stimuli
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pred_response=[]; rP_norm=[]; rO_norm=[];
        STRF_Pdc(it,1:tlag) = sub_revCOR_AM(psth_Pdc,stim_Pdc,tlag);
        
        %......................................
        % Predict response for PERIODIC stimuli
        pred_response=[]; rP_norm=[]; rO_norm=[];
%         pred_response = conv(stim_Pdc,STRF_Pdc(it,1:tlag),'full');
        pred_response = max(conv(stim_Pdc,STRF_Pdc(it,1:tlag),'full') ,0);
        
        % Normalize responses 
        rP_norm = pred_response/norm(pred_response);
        rO_norm = psth_Pdc/norm(psth_Pdc);
        
        % Correlation between predicted and observed responses 
        coinc_Pdc_Pdc(it) = max( xcorr( rP_norm, rO_norm, 0) );
        
        %.......................................
        % Predict response for IRREGULAR stimuli
        pred_response=[]; rP_norm=[]; rO_norm=[];
%         pred_response = conv(stim_Pdc,STRF_Pdc(it,1:tlag),'full');
        pred_response = max(conv(stim_Irr,STRF_Pdc(it,1:tlag),'full') ,0);
        
        % Normalize responses 
        rP_norm = pred_response/norm(pred_response);
        rO_norm = psth_Irr/norm(psth_Irr);
        
        % Correlation between predicted and observed responses 
        coinc_Pdc_Irr(it) = xcorr( rP_norm, rO_norm, 0);
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Estimate STRF from IRREGULAR stimuli
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        STRF_Irr(it,1:tlag) = sub_revCOR_AM(psth_Irr,stim_Irr,tlag);
        
        %.......................................
        % Predict response for IRREGULAR stimuli
        pred_response=[]; rP_norm=[]; rO_norm=[];
        % pred_response = conv(stim_Pdc,STRF_Pdc(it,1:tlag),'full');
        pred_response = max(conv(stim_Irr,STRF_Irr(it,1:tlag),'full') ,0);
        
        % Normalize responses
        rP_norm = pred_response/norm(pred_response);
        rO_norm = psth_Irr/norm(psth_Irr);
        
        % Correlation between predicted and observed responses
        coinc_Irr_Irr(it) = max( xcorr( rP_norm, rO_norm, 0) );
        
        %......................................
        % Predict response for PERIODIC stimuli
        pred_response=[]; rP_norm=[]; rO_norm=[];
%         pred_response = conv(stim_Pdc,STRF_Pdc(it,1:tlag),'full');
        pred_response = max(conv(stim_Pdc,STRF_Irr(it,1:tlag),'full') ,0);
        
        % Normalize responses 
        rP_norm = pred_response/norm(pred_response);
        rO_norm = psth_Pdc/norm(psth_Pdc);
        
        % Correlation between predicted and observed responses 
        coinc_Irr_Pdc(it) = xcorr( rP_norm, rO_norm, 0);
        
        
        
%                 figure;
%                 plot(rO_norm,'k')
%                 hold on
%                 plot(rP_norm,'r')
%                 title(sprintf('tlag: %i, %0.4f,  smooth win: %i',tlag,coincidence(it),smth_win))
        

    end
    
%     [c_max,i_max] = max(coinc_Pdc_Pdc);
    
%     figure;
%     plot(STRF')
    
    [c_max,i_max] = max([coinc_Pdc_Pdc coinc_Pdc_Irr coinc_Irr_Pdc coinc_Irr_Irr]);
    
    Cmax_ALL = [Cmax_ALL; c_max];
    Imax_ALL = [Imax_ALL; i_max];
    
    STRF_corr = [STRF_corr xcorr(STRF_Pdc(it,1:tlag)./norm(STRF_Pdc(it,1:tlag)),STRF_Irr(it,1:tlag)./norm(STRF_Irr(it,1:tlag)),0)];
    
    
    
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
    
    
%     Pdata_addrow = { subject session channel clu iUn  {coinc_Pdc_Pdc}  coinc_Pdc_Pdc(TLAGS == d_tlag)  c_max  i_max  range(coinc_Pdc_Pdc)  UnitData(iUn).BaseFR  mean(psth_Pdc)  VSmax  VSrt  BMFrt  };
%     Pdata = [Pdata; Pdata_addrow];
    
end %iUn


Pdata(1,:)=[];

keyboard

% How correlated are the STRFs?

% How consistent are the best tlags?
figure;
for ii=1:4
    subplot(2,2,ii);
    histogram(Cmax_ALL(:,ii),'Normalization','probability')
    xlim([0 1])
%     ylim([0 0.4])
%     set(gca,'xtick',1:11,'xticklabel',TLAGS)
end

figure;
subplot(1,2,1)
plot(Imax_ALL(:,1),Imax_ALL(:,2),'ok')
subplot(1,2,2)
plot(Imax_ALL(:,3),Imax_ALL(:,4),'ob')

    % if low correlation OR flat across tlags, check if "tuned" according
    % to VS and BMF 
    
    % check coincidence vs avg FR
    %                   vs max(VS) 
    
    % check distribution of coincidence vs PSTH binsize 



%% 

nrmRange = [Pdata.Range]./cellfun(@min,[Pdata.Corrs]);

% Plots
hf=figure;
set(hf,'Position',widescreen)

subplot(2,4,1)
histogram([Pdata.corr0],100)
xlim([0 1])
xlabel('Coinc. at tlag=100')
title('Coinc. at tlag=100')

subplot(2,4,2)
histogram([Pdata.corrBest],100)
xlim([0 1])
xlabel('Coinc. at best tlag')
title('Coinc. at best tlag')

subplot(2,4,3)
plot([0 1],[0 1],'k')
hold on
plot([Pdata.corr0],[Pdata.corrBest],'o')
plot([Pdata(nrmRange>1,:).corr0],[Pdata(nrmRange>1,:).corrBest],'ok')
plot([Pdata(Pdata.tlagBest==11,:).corr0],[Pdata(Pdata.tlagBest==11,:).corrBest],'or')
axis square
xlim([0 1])
ylim([0 1])
xlabel('Coinc. at tlag=100')
ylabel('Coinc. at best tlag')
title('Coinc100 vs Best')

subplot(2,4,4)
plot([Pdata.drivFR,],[Pdata.corr0],'o')
axis square
ylim([0 1])
xlabel('FR')
ylabel('Coinc. at tlag=100')
set(gca,'xscale','log')
title('FR vs Coinc100')


subplot(2,4,5)
histogram([Pdata.Range]./cellfun(@min,[Pdata.Corrs]),100)
hold on
histogram([Pdata(nrmRange>1,:).Range]./cellfun(@min,[Pdata(nrmRange>1,:).Corrs]),100,'FaceColor','k','FaceAlpha',1)
% xlim([0 1])
xlabel('Range, normalized to min val')
title('Range across tlags 20-500ms')

subplot(2,4,6)
histogram([Pdata.tlagBest])
set(gca,'xtick',1:length(TLAGS),'xticklabel',TLAGS)
xlabel('tlag')
title('Best tlag')

subplot(2,4,7)
plot([Pdata.maxVS,],[Pdata.corr0],'o')
axis square
xlim([0 1])
ylim([0 1])
xlabel('max VS')
ylabel('Coinc. at tlag=100')
title('max VS vs Coinc100')

subplot(2,4,8)
    yyaxis left
plot([Pdata.tlagBest]-0.1,[Pdata.maxVSrt],'o')
ylim([1 64])
set(gca,'yscale','log','ytick',AMrates,'xtick',1:length(TLAGS),'xticklabel',TLAGS)
ylabel('AMrate of max VS')
    yyaxis right
plot([Pdata.tlagBest]+0.1,[Pdata.BMFrt],'^')
axis square
ylim([1 64])
set(gca,'yscale','log','ytick',AMrates,'xtick',1:length(TLAGS),'xticklabel',TLAGS)
ylabel('BMF FR')
xlabel('Best tlag')
title('Tuning vs Best tlag')

suptitle(sprintf('smoothing win %ims',smth_win));


savedir = fullfile(fn.figs,'TRF','testing');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

savename = sprintf('Results_trPteIR_0pad_bin%i',smth_win);

print_eps_kp(hf,fullfile(savedir,savename));


end


