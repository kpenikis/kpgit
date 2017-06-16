function standardPd_FRbyUnit(subject,CRITERION)

global fn 

%~~~~~~~~~~~~~~
binsize = 25;
onlySU = 0;
%~~~~~~~~~~~~~~

fn = set_paths_directories;

% Load data
load(fullfile(fn.processed,'StandardPd_Spikes'));
InfoTable = readtable(fullfile(fn.processed,sprintf('StandardPd_StimInfo_%s',subject)));
UnitTable = readtable(fullfile(fn.processed,sprintf('UnitDataSummary_%s',subject)));


%%

% Set up figures
close all
PlotColors = [0.9 0.4 0.4; 0 0 0];
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );


% direct comparison

% FRpdc vs FRjit
hf_FRpdcFRjit(1) = figure; hold on
plot([-100 100],[-100 100],'k','LineWidth',0.5)
xlim([0 100])
ylim([0 100])
title(['FR periodic vs. FR jitter - below Crit' CRITERION])
xlabel('FR periodic')
ylabel('FR jitter')

hf_FRpdcFRjit(2) = figure; hold on
plot([-100 100],[-100 100],'k','LineWidth',0.5)
xlim([0 100])
ylim([0 100])
title(['FR periodic vs. FR jitter - above Crit' CRITERION])
xlabel('FR periodic')
ylabel('FR jitter')

% 100 ms
hf_FRprev100FRjit(1) = figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-1 1],'k','LineWidth',0.5)
xlim([-0.8 0.8])
ylim([-0.8 0.8])
title(['below ' CRITERION ' crit -- jitter FR idx vs FR idx prev 100 ms'])
xlabel('diff FR in 100 ms preceding standard pd')
ylabel('FR standard pd, jitter - periodic')

hf_FRprev100FRjit(2) = figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-1 1],'k','LineWidth',0.5)
xlim([-0.8 0.8])
ylim([-0.8 0.8])
title(['above ' CRITERION ' crit -- jitter FR idx vs FR idx prev 100 ms'])
xlabel('diff FR in 100 ms preceding standard pd')
ylabel('jitter - periodic FR idx')

% 250 ms
hf_FRprev250FRjit(1) = figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-1 1],'k','LineWidth',0.5)
xlim([-0.5 0.5])
ylim([-0.5 0.5])
title(['below ' CRITERION ' crit -- jitter FR idx vs FR idx prev 250 ms'])
xlabel('diff FR in 250 ms preceding standard pd')
ylabel('jitter - periodic FR idx')

hf_FRprev250FRjit(2) = figure; hold on
plot([-1 1],[0 0],'k','LineWidth',0.5)
plot([0 0],[-1 1],'k','LineWidth',0.5)
xlim([-0.5 0.5])
ylim([-0.5 0.5])
title(['above ' CRITERION ' crit -- jitter FR idx vs FR idx prev 250 ms'])
xlabel('diff FR in 250 ms preceding standard pd')
ylabel('jitter - periodic FR idx')

% prev rate
hf_prevPdFRjit(1) = figure; hold on
plot([0 7],[0 0],'k','LineWidth',0.5)
xlim([0 7])
ylim([-0.5 0.5])
title(['below ' CRITERION ' crit -- jitter FR idx vs Rate prev pd'])
xlabel('AM rate of previous period')
ylabel('jitter - periodic FR idx')

hf_prevPdFRjit(2) = figure; hold on
plot([0 7],[0 0],'k','LineWidth',0.5)
xlim([0 7])
ylim([-0.5 0.5])
title(['above ' CRITERION ' crit -- jitter FR idx vs Rate prev pd'])
xlabel('AM rate of previous period')
ylabel('jitter - periodic FR idx')


% Fano Factor
hf_FFpdcFFjit(1) = figure; hold on
plot([0 10],[0 10],'k','LineWidth',0.5)
xlim([0 10])
ylim([0 10])
title(['FF periodic vs. FF jitter - below Crit' CRITERION])
xlabel('FF periodic')
ylabel('FF jitter')

hf_FFpdcFFjit(2) = figure; hold on
plot([0 5],[0 5],'k','LineWidth',0.5)
xlim([0 5])
ylim([0 5])
title(['FF periodic vs. FF jitter - above Crit' CRITERION])
xlabel('FF periodic')
ylabel('FF jitter')



%% Go through data

sessions = unique(InfoTable.Session,'stable');
for isess = sessions'
    for iu = unique(InfoTable(strcmp(InfoTable.Session,isess{:}),:).cluname,'stable')'
        
        if onlySU && ~any(strcmp(InfoTable(strcmp(InfoTable.Session,isess{:}) & strcmp(InfoTable.cluname,iu{:}),:).UnitType,'SU'))
            continue
        end
        
        for iiv = unique(InfoTable(strcmp(InfoTable.Session,isess{:}) & strcmp(InfoTable.cluname,iu{:}),:).DiscrmVar,'stable')'
            
            % So far, left with one unit from one session, all data of one
            % stimulus set (at all sound params)
            
            filtInfo = InfoTable(strcmp(InfoTable.Session,isess{:}) & strcmp(InfoTable.cluname,iu{:}) & strcmp(InfoTable.DiscrmVar,iiv{:}) ,:);
            
            pars = unique([filtInfo.HP filtInfo.LP filtInfo.CenterRate filtInfo.dBSPL],'rows');
            
            for ip = 1:size(pars,1)
                
                behavs = unique(filtInfo(filtInfo.HP==pars(ip,1) & filtInfo.LP==pars(ip,2) & filtInfo.CenterRate==pars(ip,3) & filtInfo.dBSPL==pars(ip,4),:).BehavState);
                
                for ib = behavs'
                    
                    % Intermediate filtered data, just for organization
                    filtInfo2 = filtInfo(filtInfo.HP==pars(ip,1) & filtInfo.LP==pars(ip,2)...
                        & filtInfo.CenterRate==pars(ip,3) & filtInfo.dBSPL==pars(ip,4) ...
                        & strcmp(filtInfo.BehavState,ib),:);
                    
                    try
                        % Get data for periodic stimulus
                        T_pdc = filtInfo2(strcmp(filtInfo2.jitter,'0'),:);
                        
                        if isempty(T_pdc), continue, end
                        % Check that just one stimulus
                        if ~numel(unique(T_pdc.trialN))==size(T_pdc,1)
                            keyboard
                        end
                        
                        % Get unit data to filter by response quality
                        UnitData = UnitTable(ismember(UnitTable(:,1:12),T_pdc(1,[1:11 14]),'rows'),:);
                        if size(UnitData,1)>1, keyboard, end
                        
                        switch CRITERION
                            case 'VS'
                                if UnitData.VS < 0.2
                                    dpcat = 1; %below
                                else
                                    dpcat = 2; %above
                                end
                            case 'Corr'
                                if isnan(UnitData.maxCorr)
                                    dpcat = 1; %below
                                else
                                    dpcat = 2; %above
                                end
                            case {'all' 'ALL'}
                                dpcat = 1;
                        end
                        
                        
                        % Calculate response for periodic stim
                        resp_pdc = AmalgamatedSpikes(ismember(InfoTable,T_pdc,'rows'),:);
                        spbin_pdc = binspikecounts(resp_pdc,binsize);
                        FRbin_pdc = mean(spbin_pdc/binsize*1000,1);
%                         mean(reshape(sum(resp_pdc,1)/size(resp_pdc,1),[size(resp_pdc,2)/25 25]),1)*1000 ;
                        
                    catch % No periodic data for this set?
                        keyboard
                    end
                    
                    % Get unique jitter values in this stim set
                    jitters = unique(filtInfo2.jitter,'stable');
                    jitters = jitters(~strcmp(jitters,'0'));
                    
                    % Compare each jitter to the periodic condition
                    for ij = jitters'
                        
                        % Get jitter data
                        T_jit = filtInfo2(strcmp(filtInfo2.jitter,ij{:}),:);
                        
                        % Check that just one stimulus
                        if ~numel(unique(T_jit.trialN))==size(T_jit,1)
                            keyboard
                        end
                        
                        % Calculate response for jitter stim
                        resp_jit = AmalgamatedSpikes(ismember(InfoTable,T_jit,'rows'),:);
                        spbin_jit = binspikecounts(resp_jit,binsize);
                        FRbin_jit = mean(spbin_jit/binsize*1000,1);
                        
                        % Skip data points that likely come from junk clusters
                        if mean(FRbin_jit)<5 || mean(FRbin_pdc)<5
                            continue
                        end
                        
                        % Currently comparing distribution of spike counts 
                        [H,P] = ttest2(spbin_pdc,spbin_jit);
                        
                        if sum(H)>4 || any(diff(find(H),2,2)==0) || mean(P,'omitnan')<0.05
                            sigdp = 1;
                        else 
                            sigdp = 2;
                        end
                                                
                        % Plot things
                        try
                            jitFRidx = (mean(FRbin_jit)-mean(FRbin_pdc))/(mean(FRbin_jit)+mean(FRbin_pdc));
                            
                            % Fano Factor
                            figure(hf_FFpdcFFjit(dpcat))
                            plot(var(sum(spbin_jit,2))/mean(sum(spbin_jit,2)),...
                                var(sum(spbin_pdc,2))/mean(sum(spbin_pdc,2)),...
                                'o','MarkerSize',8,'LineWidth',2,...
                                'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                            % FR pdc vs. FR jit
                            figure(hf_FRpdcFRjit(dpcat))
                            plot(mean(FRbin_pdc),mean(FRbin_jit),...
                                'o','MarkerSize',8,'LineWidth',2,...
                                'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                            % FR idx histogram
%                             figure(hf_FRidxhist(dpcat))
%                             plot(mean(FRbin_pdc),mean(FRbin_jit),...
%                                 'o','MarkerSize',8,'LineWidth',2,...
%                                 'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                            % FR idx against firing history
                            figure(hf_FRprev100FRjit(dpcat))
                            plot((T_jit.prev100msFR(1)-T_pdc.prev100msFR(1))/(T_jit.prev100msFR(1)+T_pdc.prev100msFR(1)),...
                                jitFRidx,...
                                'o','MarkerSize',8,'LineWidth',2,...
                                'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                            figure(hf_FRprev250FRjit(dpcat))
                            plot((T_jit.prev250msFR(1)-T_pdc.prev250msFR(1))/(T_jit.prev250msFR(1)+T_pdc.prev250msFR(1)),...
                                jitFRidx,...
                                'o','MarkerSize',8,'LineWidth',2,...
                                'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                            % FR idx against preceding AM rate
                            rv = jitter_LUT(T_jit.CenterRate(1),ij{:});
                            figure(hf_prevPdFRjit(dpcat))
                            plot(rv(3),...
                                jitFRidx,...
                                'o','MarkerSize',8,'LineWidth',2,...
                                'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                            
                        catch
                            keyboard
                        end
                        
                    end
                end
            end
        end
    end
end


%% Save figures

savedir = fullfile(fn.standardPd,'FR');
if onlySU
    savedir = fullfile(savedir,'onlySU');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

% direct comparison
set(hf_FRpdcFRjit(1),'PaperOrientation','landscape');
print(hf_FRpdcFRjit(1),fullfile(savedir,['pdcFRjitFR_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_FRpdcFRjit(2),'PaperOrientation','landscape');
print(hf_FRpdcFRjit(2),fullfile(savedir,['pdcFRjitFR_aboveCrit' CRITERION]),'-dpdf','-bestfit');

% Fano Factor
set(hf_FFpdcFFjit(1),'PaperOrientation','landscape');
print(hf_FFpdcFFjit(1),fullfile(savedir,['pdcFFjitFF_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_FFpdcFFjit(2),'PaperOrientation','landscape');
print(hf_FFpdcFFjit(2),fullfile(savedir,['pdcFFjitFF_aboveCrit' CRITERION]),'-dpdf','-bestfit');


% previous 100 ms
set(hf_FRprev100FRjit(1),'PaperOrientation','landscape');
print(hf_FRprev100FRjit(1),fullfile(savedir,['jitFRidx_100msIdx_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_FRprev100FRjit(2),'PaperOrientation','landscape');
print(hf_FRprev100FRjit(2),fullfile(savedir,['jitFRidx_100msIdx_aboveCrit' CRITERION]),'-dpdf','-bestfit');

% previous 250 ms
set(hf_FRprev250FRjit(1),'PaperOrientation','landscape');
print(hf_FRprev250FRjit(1),fullfile(savedir,['jitFRidx_250msIdx_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_FRprev250FRjit(2),'PaperOrientation','landscape');
print(hf_FRprev250FRjit(2),fullfile(savedir,['jitFRidx_250msIdx_aboveCrit' CRITERION]),'-dpdf','-bestfit');

% previous period rate
set(hf_prevPdFRjit(1),'PaperOrientation','landscape');
print(hf_prevPdFRjit(1),fullfile(savedir,['jitFRidx_prevRate_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_prevPdFRjit(2),'PaperOrientation','landscape');
print(hf_prevPdFRjit(2),fullfile(savedir,['jitFRidx_prevRate_aboveCrit' CRITERION]),'-dpdf','-bestfit');




end


