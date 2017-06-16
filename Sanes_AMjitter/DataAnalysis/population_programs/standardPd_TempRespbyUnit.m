function standardPd_TempRespbyUnit(subject,CRITERION,single_jit)

global fn 

fn = set_paths_directories;

onlySU = 1;

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

% VSpdc vs VSjit
hf_VSpdcVSjit(1) = figure; hold on
plot([-1 1],[-1 1],'k','LineWidth',0.5)
xlim([0 0.5])
ylim([0 0.5])
title(['VS periodic vs. VS jitter - below Crit' CRITERION])
xlabel('VS periodic')
ylabel('VS jitter')

hf_VSpdcVSjit(2) = figure; hold on
plot([-1 1],[-1 1],'k','LineWidth',0.5)
xlim([0 0.5])
ylim([0 0.5])
title(['VS periodic vs. VS jitter - above Crit' CRITERION])
xlabel('VS periodic')
ylabel('VS jitter')

% Rayleigh Stat pdc vs Rayleigh Stat jitter
hf_RaySpdc_RaySjit(1) = figure; hold on
plot([-1 100],[-1 100],'k','LineWidth',0.5)
xlim([0 25])
ylim([0 25])
title(['Rayleigh Statistic periodic vs. Rayleigh Statistic jitter - below Crit' CRITERION])
xlabel('Rayleigh Statistic periodic')
ylabel('Rayleigh Statistic jitter')

hf_RaySpdc_RaySjit(2) = figure; hold on
plot([-1 100],[-1 100],'k','LineWidth',0.5)
xlim([0 25])
ylim([0 25])
title(['Rayleigh Statistic periodic vs. Rayleigh Statistic jitter - above Crit' CRITERION])
xlabel('Rayleigh Statistic periodic')
ylabel('Rayleigh Statistic jitter')

% Response Strength periodic vs Response Strength jitter
hf_RespSpdc_RespSjit(1) = figure; hold on
plot([-1 5],[-1 5],'k','LineWidth',0.5)
xlim([0 5])
ylim([0 5])
title(['Resp Strength (VS * avg NSpks) periodic vs. Resp Strength jitter - below Crit' CRITERION])
xlabel('Response Strength periodic')
ylabel('Response Strength jitter')

hf_RespSpdc_RespSjit(2) = figure; hold on
plot([-1 5],[-1 5],'k','LineWidth',0.5)
xlim([0 5])
ylim([0 5])
title(['Resp Strength (VS * avg NSpks) periodic vs. Resp Strength jitter - above Crit' CRITERION])
xlabel('Response Strength periodic')
ylabel('Response Strength jitter')


% VS difference index regressed against various things

% 100 ms
hf_VSidx_100msidx(1) = figure; hold on
plot([0 0],[-1 1],'k','LineWidth',0.5)
plot([-1 1],[0 0],'k','LineWidth',0.5)
xlim([-1 1])
ylim([-1 1])
title(['VS index as a fct of normalized change prior 100 ms FR - below Crit' CRITERION])
xlabel('previous 100 ms FR idx')
ylabel('jitter-periodic VS idx')

hf_VSidx_100msidx(2) = figure; hold on
plot([0 0],[-1 1],'k','LineWidth',0.5)
plot([-1 1],[0 0],'k','LineWidth',0.5)
xlim([-1 1])
ylim([-1 1])
title(['VS index as a fct of normalized change prior 100 ms FR - above Crit' CRITERION])
xlabel('previous 100 ms FR idx')
ylabel('jitter-periodic VS idx')

hf_VSidx_100ms(1) = figure; hold on
plot([-1 50],[0 0],'k','LineWidth',0.5)
xlim([0 50])
ylim([-1 1])
title(['VS index as a fct of prior 100 ms FR - below Crit' CRITERION])
xlabel('previous 100 ms FR (JITTER)')
ylabel('jitter-periodic VS idx')

hf_VSidx_100ms(2) = figure; hold on
plot([-1 50],[0 0],'k','LineWidth',0.5)
xlim([0 50])
ylim([-1 1])
title(['VS index as a fct of prior 100 ms FR - above Crit' CRITERION])
xlabel('previous 100 ms FR (JITTER)')
ylabel('jitter-periodic VS idx')


% 250 ms
hf_VSidx_250msidx(1) = figure; hold on
plot([0 0],[-1 1],'k','LineWidth',0.5)
plot([-1 1],[0 0],'k','LineWidth',0.5)
xlim([-1 1])
ylim([-1 1])
title(['VS index as a fct of normalized change prior 250 ms FR - below Crit' CRITERION])
xlabel('previous 250 ms FR idx')
ylabel('jitter-periodic VS idx')

hf_VSidx_250msidx(2) = figure; hold on
plot([0 0],[-1 1],'k','LineWidth',0.5)
plot([-1 1],[0 0],'k','LineWidth',0.5)
xlim([-1 1])
ylim([-1 1])
title(['VS index as a fct of normalized change prior 250 ms FR - above Crit' CRITERION])
xlabel('previous 250 ms FR idx')
ylabel('jitter-periodic VS idx')

hf_VSidx_250ms(1) = figure; hold on
plot([-1 50],[0 0],'k','LineWidth',0.5)
xlim([0 50])
ylim([-1 1])
title(['VS index as a fct of prior 250 ms FR - below Crit' CRITERION])
xlabel('previous 250 ms FR (JITTER)')
ylabel('jitter-periodic VS idx')

hf_VSidx_250ms(2) = figure; hold on
plot([-1 50],[0 0],'k','LineWidth',0.5)
xlim([0 50])
ylim([-1 1])
title(['VS index as a fct of prior 250 ms FR - above Crit' CRITERION])
xlabel('previous 250 ms FR (JITTER)')
ylabel('jitter-periodic VS idx')


% Previous period rate
hf_VSidx_prevPd(1) = figure; hold on
plot([0 7],[0 0],'k','LineWidth',0.5)
xlim([0 7])
ylim([-1 1])
title(['VS index as a fct of AM rate of previous period - below Crit' CRITERION])
xlabel('AM rate of previous period')
ylabel('jitter-periodic VS idx')

hf_VSidx_prevPd(2) = figure; hold on
plot([0 7],[0 0],'k','LineWidth',0.5)
xlim([0 7])
ylim([-1 1])
title(['VS index as a fct of AM rate of previous period - above Crit' CRITERION])
xlabel('AM rate of previous period')
ylabel('jitter-periodic VS idx')


% Go through data

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
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Now we have data from one unit, one session, and one
                    % set of acoustic parameters (everything but jitters).
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    try
                        % Get data for periodic stimulus
                        T_pdc = filtInfo2(strcmp(filtInfo2.jitter,'0'),:);
                        
                        if isempty(T_pdc), continue, end
                        % Check that just one stimulus
                        if ~numel(unique(T_pdc.trialN))==size(T_pdc,1)
                            keyboard
                        end
                        
                        % Calculate response for periodic stim
                        resp_pdc = AmalgamatedSpikes(ismember(InfoTable,T_pdc,'rows'),:);
                        
                        [VS_pdc,RayS_pdc] = vectorstrength(find(reshape(resp_pdc,[1 size(resp_pdc,1)*size(resp_pdc,2)])), size(resp_pdc,2), size(resp_pdc,1));
                        RespS_pdc = VS_pdc * sum(sum(resp_pdc))/size(resp_pdc,1);
%                         meanPhase_pdc = calc_mean_phase(resp_pdc);
                        
                        spbin_pdc = binspikecounts(resp_pdc,10);
                        FRbin_pdc = mean(reshape(sum(resp_pdc,1)/size(resp_pdc,1),[size(resp_pdc,2)/25 25]),1)*1000 ;
                        
                    catch % No periodic data for this set?
                        keyboard
                    end
                    
                    %% Categorize by temporal response or not 
                    % Get unit data to filter by response quality
                    UnitData = UnitTable(ismember(UnitTable(:,1:12),T_pdc(1,[1:11 14]),'rows'),:);
                    if size(UnitData,1)>1, keyboard, end
                    
                    use_crit = UnitData.VS;
%                     use_crit = VS_pdc;
                    
                    switch CRITERION
                        case 'VS'
                            if use_crit < 0.2
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
                        case 'ALL'
                            dpcat = 1;
                    end
                    
                    %% Get jitters and plot
                    % Get unique jitter values in this stim set
                    jitters = unique(filtInfo2.jitter,'stable');
                    jitters = jitters(~strcmp(jitters,'0'));
                    
                    % Compare each jitter to the periodic condition
                    for ij = jitters'
                        
                        if ~strcmpi(single_jit,'ALL') && ~strcmp(ij{:},single_jit)
                            continue
                        end
                        
                        % Get jitter data
                        T_jit = filtInfo2(strcmp(filtInfo2.jitter,ij{:}),:);
                        % Check that just one stimulus
                        if ~(numel(unique(T_jit.trialN))==size(T_jit,1))
                            keyboard
                        end
                        
                        % Get spikes for jitter stim
                        resp_jit = AmalgamatedSpikes(ismember(InfoTable,T_jit,'rows'),:);
                        
                        % Calculate response for jitter stim
                        [VS_jit,RayS_jit] = vectorstrength(find(reshape(resp_jit,[1 size(resp_jit,1)*size(resp_pdc,2)])), size(resp_jit,2), size(resp_jit,1));
                        RespS_jit = VS_jit * sum(sum(resp_jit))/size(resp_jit,1);
                        
                        spbin_jit = binspikecounts(resp_jit,10);
                        FRbin_jit = mean(reshape(sum(resp_jit,1)/size(resp_jit,1),[size(resp_jit,2)/25 25]),1)*1000 ;
                        
                        % Skip data points that likely come from junk clusters
                        if mean(FRbin_jit)<5 || mean(FRbin_pdc)<5
                            continue
                        end
                        
%                         % Currently comparing distribution of spike counts 
%                         [H,P] = ttest2(spbin_pdc,spbin_jit);
%                         
%                         if sum(H)>4 || any(diff(find(H),2,2)==0) || mean(P,'omitnan')<0.05
%                             sigdp = 1;
%                         else 
%                             sigdp = 2;
%                         end
                        
                        % Plot things
                        try
                            
                        % Periodic to Jitter Comparisons
                        figure(hf_VSpdcVSjit(dpcat))
                        plot(VS_pdc,...
                            VS_jit,...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        figure(hf_RaySpdc_RaySjit(dpcat))
                        plot(RayS_pdc,...
                            RayS_jit,...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        figure(hf_RespSpdc_RespSjit(dpcat))
                        plot(RespS_pdc,...
                            RespS_jit,...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        % Regressed variables
                        figure(hf_VSidx_100ms(dpcat))
                        plot(T_jit.prev100msFR(1),...
                            (VS_jit-VS_pdc)/(VS_jit+VS_pdc),...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        figure(hf_VSidx_100msidx(dpcat))
                        plot((T_jit.prev100msFR(1)-T_pdc.prev100msFR(1))/(T_jit.prev250msFR(1)+T_pdc.prev250msFR(1)),...
                            (VS_jit-VS_pdc)/(VS_jit+VS_pdc),...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        figure(hf_VSidx_250ms(dpcat))
                        plot(T_jit.prev250msFR(1),...
                            (VS_jit-VS_pdc)/(VS_jit+VS_pdc),...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                                                
                        figure(hf_VSidx_250msidx(dpcat))
                        plot((T_jit.prev250msFR(1)-T_pdc.prev250msFR(1))/(T_jit.prev250msFR(1)+T_pdc.prev250msFR(1)),...
                            (VS_jit-VS_pdc)/(VS_jit+VS_pdc),...
                            'o','MarkerSize',8,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));%PlotColors(sigdp,:));
                        
                        figure(hf_VSidx_prevPd(dpcat))
                        rv = jitter_LUT(T_jit.CenterRate(1),ij{:});
                        plot(rv(3),...
                            (VS_jit-VS_pdc)/(VS_jit+VS_pdc),...
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

savedir = fullfile(fn.standardPd,'Temporal');
if onlySU
    savedir = fullfile(savedir,'onlySU');
end
if ~strcmpi(single_jit,'ALL')
    savedir = fullfile(savedir,single_jit);
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

if strcmpi(CRITERION,'ALL')
    crit = {CRITERION};
else
    crit = {['below' CRITERION] ['above' CRITERION]};
end

for ic = 1:numel(crit)
    
    % Vector Strength
    set(hf_VSpdcVSjit(ic),'PaperOrientation','landscape');
    print(hf_VSpdcVSjit(ic),fullfile(savedir,['VS_comparison_' crit{ic}]),'-dpdf','-bestfit');
    
    % Rayleigh Statistic
    set(hf_RaySpdc_RaySjit(ic),'PaperOrientation','landscape');
    print(hf_RaySpdc_RaySjit(ic),fullfile(savedir,['RayleighStat_comparison_' crit{ic}]),'-dpdf','-bestfit');
    
    % Response Strength (VS*NSpk)
    set(hf_RespSpdc_RespSjit(ic),'PaperOrientation','landscape');
    print(hf_RespSpdc_RespSjit(ic),fullfile(savedir,['ResponseStrength_comparison' crit{ic}]),'-dpdf','-bestfit');
    
    
    % jit-pdc VS index
    set(hf_VSidx_100ms(ic),'PaperOrientation','landscape');
    print(hf_VSidx_100ms(ic),fullfile(savedir,['VSidx_100msJit' crit{ic}]),'-dpdf','-bestfit');
    
    set(hf_VSidx_100msidx(ic),'PaperOrientation','landscape');
    print(hf_VSidx_100msidx(ic),fullfile(savedir,['VSidx_100msidx' crit{ic}]),'-dpdf','-bestfit');
    
    set(hf_VSidx_250ms(ic),'PaperOrientation','landscape');
    print(hf_VSidx_250ms(ic),fullfile(savedir,['VSidx_250msJit' crit{ic}]),'-dpdf','-bestfit');
    
    set(hf_VSidx_250msidx(ic),'PaperOrientation','landscape');
    print(hf_VSidx_250msidx(ic),fullfile(savedir,['VSidx_250msidx' crit{ic}]),'-dpdf','-bestfit');
    
    set(hf_VSidx_prevPd(ic),'PaperOrientation','landscape');
    print(hf_VSidx_prevPd(ic),fullfile(savedir,['VSidx_prevPd' crit{ic}]),'-dpdf','-bestfit');
    
end



end


