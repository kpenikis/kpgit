function TRF_PreObs_AMVS(RERUN)
%
%  TRF_PreObs_AMVS
%    version with splitting trials and bootstrapping 
%    fit model with AM data, test on speech 
%
%  KP, 2019-11
%

%Xadd durations of VS stimuli
%Xmatch AM and VS seconds of data? NO DIFFERENCE IN RESULT
%Xmatch amplitude of stimuli across sessions
% Range of raw sound recorded is different in AM and VS of Apr11
%  AM: rms=0.5779, min=-3.1794, max=3.2113
%  VS: rms=0.0104, min=-0.0997, max=0.0987
% Check RPVDS; think about stimulus differences; listen
% FOR NOW, normalize both to be in same range
% NO DIFFERENCE WITH RAW OR NORMALIZED STIMULI
%Xcompare predictions to within context predictions



close all
global fn spkshift smth_win exclOnset AM_durs VS_durs nTrGrp

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
smth_win    = 5;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% only set up for one val rn
TLAGS       = 500;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
exclOnset   = 0; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PlotEx      = 0;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nBS         = 50;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nTrGrp      = 10;

if nargin<1
    RERUN = 1;
end


%% Units setup

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UData_AM = q.UnitData;
UInfo_AM = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UData_AM([UData_AM.IntTime_spk]>0).IntTime_spk]);
%-------

q = load(fullfile(fn.processed,'UnitsVS'));
UData_VS = q.UnitData;
UInfo_VS = q.UnitInfo;
clear q

% Load AM-VS LUT
load(fullfile(fn.processed,'UnMatchedLUT'));

Uindices_AMVS = nan(size(MatchedUnits,1),2);
for iu = 1:size(MatchedUnits,1)
    iU_AM = []; iU_VS = []; 
    iU_AM = find( strcmp(MatchedUnits.AMsessions{iu},{UData_AM.Session}) & MatchedUnits.AMclu(iu)==[UData_AM.Clu] );
    iU_VS = find( strcmp(MatchedUnits.VSsessions{iu},{UData_VS.Session}) & MatchedUnits.VSclu(iu)==[UData_VS.Clu] );
    if numel(iU_AM)==1 && numel(iU_VS)==1
        Uindices_AMVS(iu,:) = [iU_AM iU_VS];
    else
        keyboard
    end
end


%% Data setup

AMrates   = [2 4 8 16 32];
AM_durs = [1500   1000   1000   1000   1000   1000   1937   1937   1000];
VS_durs = [1670   1406   2412   5928   2554   2556];

rng('shuffle')

% Preallocate
Results = table;
Results.iiUn          = nan;
Results.smthWin       = nan;
Results.tlag          = nan;
Results.C_AM_self     = nan;
Results.C_AM_self_std = nan;
Results.C_AM_VS       = nan;
Results.C_AM_VS_std   = nan;
Results.C_VS_self     = nan;
Results.C_VS_self_std = nan;
Results.C_VS_AM       = nan;
Results.C_VS_AM_std   = nan;
Results.BaseFR        = nan;
Results.BMFrt         = nan;
Results.UnType        = {''};


%% Figure settings

savedir = fullfile(fn.figs,'TRF_AMVS','Feb2020');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallrect    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

dotSize = 30;


%%

if RERUN

% %     STIM_AM_ALL = cell(7,1);
% %     STIM_VS_ALL = cell(7,1);
    
for iiUn = 1:size(Uindices_AMVS,1)
    
    iUnAM = Uindices_AMVS(iiUn,1);
    iUnVS = Uindices_AMVS(iiUn,2);
    
    UnType = '';
    if UInfo_AM(iUnAM,:).TroughPeak>0.43 &&  UInfo_VS(iUnVS,:).TroughPeak>0.43
        UnType = 'RS';
    elseif UInfo_AM(iUnAM,:).TroughPeak<=0.43 &&  UInfo_VS(iUnVS,:).TroughPeak<=0.43
        UnType = 'NS';
    else
        fprintf('\nskipping iiUn %i \n\n',iiUn)
        skipUns = [skipUns; iiUn];
        continue
    end
    
    
    % Get PSTH and Stim data
    [PSTH_AM,STIM_AM] = getdata4TRF(UData_AM,iUnAM);
    [PSTH_VS,STIM_VS] = getdata4TRF(UData_VS,iUnVS);
    
% % %     % Save averages for fft later
% % %     for ist = 1:numel(STIM_AM)
% % %         STIM_AM_ALL{ist} = [STIM_AM_ALL{ist}; mean(STIM_AM{ist},1)];
% % %     end
% % %     for ist = 1:numel(STIM_VS)
% % %         STIM_VS_ALL{ist} = [STIM_VS_ALL{ist}; mean(STIM_VS{ist},1)];
% % %     end
    
    fprintf(' running model\n')
     
    
    %%  ------  STRF analysis  ------
    
    if PlotEx
        hfex=figure;
        set(hfex,'Position',fullscreen,'NextPlot','add')
    end
    
    nTrAll = cellfun(@(x) size(x,1), [PSTH_AM; PSTH_VS]);
    
    STRF_AM       = nan(numel(TLAGS),max(TLAGS),nBS);
    STRF_VS       = nan(numel(TLAGS),max(TLAGS),nBS);
    coinc_AM_self = nan(numel(TLAGS),nBS);
    coinc_AM_VS   = nan(numel(TLAGS),nBS);
    coinc_VS_self = nan(numel(TLAGS),nBS);
    coinc_VS_AM   = nan(numel(TLAGS),nBS);
    
    savetimewin = 5400 + (1:6000);
    rP_diff_SAVE = [];
    rP_self_SAVE = [];
    rObs_SAVE    = [];
    
    for iBS = 1:nBS
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Concatenate all AM data
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        psth_AM1 = []; stim_AM1 = [];
        psth_AM2 = []; stim_AM2 = [];
        
        for st = find(~cellfun(@isempty,PSTH_AM))'
            
            if size(STIM_AM{st},1)<(nTrGrp*2)
                continue
            end
            
            allTrs = randperm(size(STIM_AM{st},1));
            
            % Average trials and concatenate stim
            psth_AM1  = [psth_AM1 mean(PSTH_AM{st}(allTrs(1:nTrGrp),:),1) ];
            stim_AM1  = [stim_AM1 mean(STIM_AM{st}(allTrs(1:nTrGrp),:),1) ];
            
            psth_AM2  = [psth_AM2 mean(PSTH_AM{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
            stim_AM2  = [stim_AM2 mean(STIM_AM{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
            
        end
        if length(psth_AM1)<5000
            keyboard
        end
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Concatenate SPEECH data
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        psth_VS1 = []; stim_VS1 = [];
        psth_VS2 = []; stim_VS2 = [];
        stim_plot = [];
        
        for st = find(~cellfun(@isempty,PSTH_VS))'
            
            if size(STIM_VS{st},1)<(nTrGrp*2)
                continue
            end
            
            allTrs = randperm(size(STIM_VS{st},1));
            
            % Average trials and concatenate stim
            psth_VS1  = [psth_VS1 mean(PSTH_VS{st}(allTrs(1:nTrGrp),:),1) ];
            stim_VS1  = [stim_VS1 mean(STIM_VS{st}(allTrs(1:nTrGrp),:),1) ];
            
            psth_VS2  = [psth_VS2 mean(PSTH_VS{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
            stim_VS2  = [stim_VS2 mean(STIM_VS{st}(allTrs(nTrGrp+(1:nTrGrp)),:),1) ];
            
            stim_plot = [stim_plot mean(STIM_VS{st}(allTrs,:),1) ];
        end
        
%         figure;
%         subplot(4,1,1);
%         hold on
%         plot(stim_plot(savetimewin),'k','LineWidth',2)
%         xlim([0 size(savetimewin,2)])
%         print_eps_kp(gcf,fullfile(savedir,'ExPSTH_stim'))
        
        if length(psth_VS1)<2000
            aaa=-245;
        end
        
        
        % Normalize stimulus traces to [0 1]
        minAM    = min([stim_AM1 stim_AM2]);
        maxAM    = max([stim_AM1-minAM stim_AM2-minAM]);
        stim_AM1 = (stim_AM1-minAM)./maxAM;
        stim_AM2 = (stim_AM2-minAM)./maxAM;
        
        minVS    = min([stim_VS1 stim_VS2]);
        maxVS    = max([stim_VS1-minVS stim_VS2-minVS]);
        stim_VS1 = (stim_VS1-minVS)./maxVS;
        stim_VS2 = (stim_VS2-minVS)./maxVS;
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %  Fit STRFs and compare responses
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % For each tlag, calculate STRF, predict response, and compare to observed
        
        for it = 1:numel(TLAGS)
            
            tlag  = TLAGS(it);
            
            %......................................
            %                TRAIN
            %......................................
            
            % Estimate STRF from AM stimuli
            STRF_AM(it,1:tlag,iBS) = sub_revCOR_AM(psth_AM1,stim_AM1,tlag);
            
            
            % Estimate STRF from VS stimuli
            STRF_VS(it,1:tlag,iBS) = sub_revCOR_AM(psth_VS1,stim_VS1,tlag);
            
            
            %......................................
            %                 TEST
            %......................................
            
            %~~~~~~~~~~~~~~
            %  Predict AM
            %~~~~~~~~~~~~~~
            %- - - - - - -  from VS  - - - - - - -
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_AM2,STRF_VS(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_AM2/norm(psth_AM2);
            
            % Correlation between predicted and observed responses
            coinc_AM_VS(it,iBS) = max( xcorr( rP_norm, rO_norm, 0) );
            
            if PlotEx
%                 figure(hfex); clf
%                 subplot(2,2,1);
%                 hold on
%                 plot(rO_norm,'k-')
%                 plot(rP_norm,'r-','LineWidth',2)
%                 xlim([0 length(rO_norm)])
%                 title(sprintf('Test AM, train VS, coinc=%0.3f',coinc_AM_VS(it,iBS)))
            end
            
            %- - - - - - - from self  - - - - - - -
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_AM2,STRF_AM(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_AM2/norm(psth_AM2);
            
            % Correlation between predicted and observed responses
            coinc_AM_self(it,iBS) = max( xcorr( rP_norm, rO_norm, 0) );
            
            if PlotEx
%                 subplot(2,2,3);
%                 hold on
%                 plot(rO_norm,'k-')
%                 plot(rP_norm,'b-','LineWidth',2)
%                 xlim([0 length(rO_norm)])
%                 title(sprintf('Test AM, train self, coinc=%0.3f',coinc_AM_self(it,iBS)))
            end
            
            
            %~~~~~~~~~~~~~~
            %  Predict VS
            %~~~~~~~~~~~~~~
            %- - - - - - -  from AM  - - - - - - -
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_VS2,STRF_AM(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_VS2/norm(psth_VS2);
            
            % Save for plot later
            rP_diff_SAVE = [rP_diff_SAVE; rP_norm(savetimewin)];
            
            % Correlation between predicted and observed responses
            coinc_VS_AM(it,iBS) = max( xcorr( rP_norm, rO_norm, 0) );
            
            if PlotEx
                subplot(2,1,1);
                hold on
                plot(rO_norm,'k-')
                plot(rP_norm,'r-','LineWidth',2)
                xlim([0 length(rO_norm)])
                title(sprintf('Test VS, train AM, coinc=%0.3f',coinc_VS_AM(it,iBS)))
            end
            
            %- - - - - - - from self  - - - - - - -
            pred_response=[]; rP_norm=[]; rO_norm=[];
            pred_response = max(conv(stim_VS2,STRF_VS(it,1:tlag,iBS),'full') ,0);
            
            % Normalize responses
            rP_norm = pred_response/norm(pred_response);
            rO_norm = psth_VS2/norm(psth_VS2);
            
            % Save for plot later
            rP_self_SAVE = [rP_self_SAVE; rP_norm(savetimewin)];
            rObs_SAVE = [rObs_SAVE; rO_norm(savetimewin)];
            
            % Correlation between predicted and observed responses
            coinc_VS_self(it,iBS) = max( xcorr( rP_norm, rO_norm, 0) );
            
            if PlotEx
                subplot(2,1,2);
                hold on
                plot(rO_norm,'k-')
                plot(rP_norm,'b-','LineWidth',2)
                xlim([0 length(rO_norm)])
                title(sprintf('Test VS, train self, coinc=%0.3f',coinc_VS_self(it,iBS)))
            end
            
            %             print_eps_kp(hfex,fullfile(savedir,['ExamplePredObs_iUn' num2str(iUnAM)]))
            
            
        end %tlags
        
%         plot(STRF_AM(it,1:tlag,iBS),'m')
%         plot(STRF_VS(it,1:tlag,iBS),'b')
        
    end %iBS
    
    figure;
    
    subplot(2,1,1);
    hold on
    plot(mean(rObs_SAVE,1),'k-')
    plot(mean(rP_diff_SAVE,1),'r-','LineWidth',2)
    xlim([0 length(rObs_SAVE)])
    title(sprintf('Test VS, train AM, coinc=%0.3f',mean(coinc_VS_AM(it,:))))
    
    subplot(2,1,2);
    hold on
    plot(mean(rObs_SAVE,1),'k-')
    plot(mean(rP_self_SAVE,1),'b-','LineWidth',2)
    xlim([0 length(rObs_SAVE)])
    title(sprintf('Test VS, train self, coinc=%0.3f',mean(coinc_VS_self(it,:))))
    
    keyboard
    print_eps_kp(gcf,fullfile(savedir,sprintf('ExPSTH_%i',iiUn)))
    
    
%     figure;
%     histogram(coinc_VS_AM)
%     hold on
%     histogram(coinc_AM_self)
%     histogram(coinc_AM_VS)
%     histogram(coinc_VS_self)
%     legend({'VS from AM' 'AM self' 'AM from VS' 'VS self'},'Location','best')
%     title([MatchedUnits.AMsessions{iiUn} ' ' num2str(MatchedUnits{iiUn,2})...
%         ', ' MatchedUnits.VSsessions{iiUn} ' ' num2str(MatchedUnits{iiUn,4})])
%     xlim([0 0.6])
%     
%     savename = sprintf('CoincVals_%i_long',iiUn);
%     print_eps_kp(gcf,fullfile(savedir,savename))
%     
%     
    
    % if more than 50% nans (bc few spikes) just skip
    if mode(sum(isnan(coinc_VS_AM),2)/size(coinc_VS_AM,2)) > 0.5
        keyboard
        coinc_VS_AM = nan(size(coinc_VS_AM));
    end
    
    
    %% Save data to table
    
    BMFrt = AMrates(UData_AM(iUnAM).iBMF_FR);
    if isempty(BMFrt)
        BMFrt=nan;
    end
    
    % Add to Results table
    Results_addrow = { iiUn smth_win tlag  ...
        mean(coinc_AM_self,'omitnan') std(coinc_AM_self,'omitnan')...
        mean(coinc_AM_VS,'omitnan')   std(coinc_AM_VS,'omitnan')...
        mean(coinc_VS_self,'omitnan') std(coinc_VS_self,'omitnan')...
        mean(coinc_VS_AM,'omitnan')   std(coinc_VS_AM,'omitnan')...
        UData_AM(iUnAM).BaseFR   BMFrt  UnType  };
    Results = [Results; Results_addrow];
    
    
end %iUn

Results(1,:)=[];

savename = sprintf('Res_TRF_%i_AMVS_On-%i',TLAGS,sum(~exclOnset));
% Save Results
save(fullfile(savedir,savename),'Results','-v7.3')

% And write to Matched Units excel file
% writetable(Results,fullfile(fn.processed,'UnMatchedLUT.xlsx'),...
%     'Sheet',4,'Range','A1:K30')


else
    
    % Load Results
    savename = sprintf('Res_TRF_%i_AMVS_On-%i',TLAGS,sum(~exclOnset));
    q=load(fullfile(savedir,savename));
    Results = q.Results;
    
    
end


%% Plot results

axmax = 1; 
iRS = strcmp(Results.UnType,'RS');
iNS = strcmp(Results.UnType,'NS');

BMs = sum(MatchedUnits{Results.iiUn,end-2:end},2)<4;


%.......................................
%   New figure
%.......................................

hf=figure;
set(hf,'Position',tallrect,'NextPlot','add')

%~~~~~~~~~~
%   AM
%~~~~~~~~~~
hs(1)=subplot(2,2,1);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

% plot([Results.C_AM_self Results.C_VS_self]',[Results.C_AM_VS Results.C_VS_AM]',...
%     '-k','LineWidth',2)
plot(Results(iRS,:).C_AM_self,Results(iRS,:).C_AM_VS,'.r','MarkerSize',dotSize)
plot(Results(iNS,:).C_AM_self,Results(iNS,:).C_AM_VS,'.b','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title(sprintf('All putative matches, tlag=%i (N=%i)',TLAGS,size(Results,1)))
xlabel('Same Context (AM)')
ylabel('Opposite Context')

% Stats
C_Sam = Results.C_AM_self;
C_Opp = Results.C_AM_VS;
ig = ~isnan(C_Sam) & ~isnan(C_Opp);

p_wsr = signrank(C_Opp(ig),C_Sam(ig));
avg_diff = mean(C_Opp(ig)-C_Sam(ig),'omitnan');

[r_c,p_c] = corr(C_Opp(ig),C_Sam(ig));

text(0.02,0.9,sprintf('signrank p=%0.2e\navg diff C = %0.4f\nPearson r=%0.2f, p=%0.2e',p_wsr,avg_diff,r_c,p_c))


%~~~~~~~~~~
%   VS
%~~~~~~~~~~
hs(2)=subplot(2,2,2);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

% plot([Results(iRS,:).C_AM_self Results.C_AM_VS]',[Results.C_VS_self Results.C_VS_AM]',...
%     '.-k','LineWidth',2,'MarkerSize',dotSize,'MarkerFaceColor','k')
plot(Results(iRS,:).C_VS_self,Results(iRS,:).C_VS_AM,'.r','MarkerSize',dotSize)
plot(Results(iNS,:).C_VS_self,Results(iNS,:).C_VS_AM,'.b','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title('')
xlabel('Same Context (VS)')
ylabel('Opposite Context')

% Stats
C_Sam = Results.C_VS_self;
C_Opp = Results.C_VS_AM;
ig = ~isnan(C_Sam) & ~isnan(C_Opp);

p_wsr = signrank(C_Opp(ig),C_Sam(ig));
avg_diff = mean(C_Opp(ig)-C_Sam(ig),'omitnan');

[r_c,p_c] = corr(C_Opp(ig),C_Sam(ig));

text(0.02,0.9,sprintf('signrank p=%0.2e\navg diff C = %0.4f\nPearson r=%0.2f, p=%0.2e',p_wsr,avg_diff,r_c,p_c))



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Just BEST matches
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% AM
hs(3)=subplot(2,2,3);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

% plot([Results(BMs,:).C_AM_self Results(BMs,:).C_VS_self]',[Results(BMs,:).C_AM_VS Results(BMs,:).C_VS_AM]',...
%     '-k','LineWidth',2)
plot(Results(BMs&iRS,:).C_AM_self,Results(BMs&iRS,:).C_AM_VS,'.r','MarkerSize',dotSize)
plot(Results(BMs&iNS,:).C_AM_self,Results(BMs&iNS,:).C_AM_VS,'.b','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title(sprintf('Best matches, N=%i',sum(BMs)))
xlabel('Same Context (AM)')
ylabel('Opposite Context')


% VS
hs(4)=subplot(2,2,4);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

% plot([Results(BMs,:).C_AM_self Results(BMs,:).C_VS_self]',[Results(BMs,:).C_AM_VS Results(BMs,:).C_VS_AM]',...
%     '-k','LineWidth',2)
plot(Results(BMs&iRS,:).C_VS_self,Results(BMs&iRS,:).C_VS_AM,'.r','MarkerSize',dotSize)
plot(Results(BMs&iNS,:).C_VS_self,Results(BMs&iNS,:).C_VS_AM,'.b','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title(sprintf('Best matches, N=%i',sum(BMs)))
xlabel('Same Context (VS)')
ylabel('Opposite Context')


% Save figure
print_eps_kp(hf,fullfile(savedir,savename))



keyboard




%.......................................
% Original figure
%.......................................

hf=figure;
set(hf,'Position',tallrect,'NextPlot','add')

% Scatter plot
hs(1)=subplot(2,2,1);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

plot([Results.C_AM_self Results.C_VS_self]',[Results.C_AM_VS Results.C_VS_AM]',...
    '-k','LineWidth',2)
plot(Results.C_AM_self,Results.C_AM_VS,'.m','MarkerSize',dotSize,'MarkerEdgeColor','m')
plot(Results.C_VS_self,Results.C_VS_AM,'.g','MarkerSize',dotSize,'MarkerFaceColor','g')

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title(sprintf('All putative matches, tlag=%i (N=%i)',tlag,size(Results,1)))
xlabel('Same Context Prediction')
ylabel('Opposite Context Prediction')


% AM vs VS
hs(2)=subplot(2,2,2);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

plot([Results.C_AM_self Results.C_AM_VS]',[Results.C_VS_self Results.C_VS_AM]',...
    '.-k','LineWidth',2,'MarkerSize',dotSize,'MarkerFaceColor','k')
plot(Results.C_AM_VS,Results.C_VS_AM,'b.','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title('')
xlabel('AM Context Prediction')
ylabel('VS Context Prediction')




% Now, just best cells and matches

BMs = sum(MatchedUnits{:,end-2:end},2)==0;


% Scatter plot
hs(3)=subplot(2,2,3);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

plot([Results(BMs,:).C_AM_self Results(BMs,:).C_VS_self]',[Results(BMs,:).C_AM_VS Results(BMs,:).C_VS_AM]',...
    '-k','LineWidth',2)
plot(Results(BMs,:).C_AM_self,Results(BMs,:).C_AM_VS,'.m','MarkerSize',dotSize,'MarkerEdgeColor','m')
plot(Results(BMs,:).C_VS_self,Results(BMs,:).C_VS_AM,'.g','MarkerSize',dotSize,'MarkerFaceColor','g')

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title(sprintf('Best matches, N=%i',sum(BMs)))
xlabel('Same Context Prediction')
ylabel('Opposite Context Prediction')


% AM vs VS
hs(4)=subplot(2,2,4);
hold on
plot([0 axmax],[0 axmax],'Color',[0.7 0.7 0.7])
axis square; box on

plot([Results(BMs,:).C_AM_self Results(BMs,:).C_AM_VS]',[Results(BMs,:).C_VS_self Results(BMs,:).C_VS_AM]',...
    '.-k','LineWidth',2,'MarkerSize',dotSize,'MarkerFaceColor','k')
plot(Results(BMs,:).C_AM_VS,Results(BMs,:).C_VS_AM,'b.','MarkerSize',dotSize)

xlim([0 axmax])
ylim([0 axmax])
set(gca,'Color','none')
title('')
xlabel('AM Context Prediction')
ylabel('VS Context Prediction')


% Save figure
print_eps_kp(hf,fullfile(savedir,savename))




end


