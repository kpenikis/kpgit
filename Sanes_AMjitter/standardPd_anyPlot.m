function standardPd_anyPlot(subject,plottype,var1,var2,categorization)

global fn Figs


% Plot periodic data on own fig
% remove formatting from title
% try FR or peak, now that restricting n trials


%~~~~~~~~~~~~~~
binsize = 25;
onlySU = 1;
trThresh = 18;
%~~~~~~~~~~~~~~

fn = set_paths_directories;

% Load data
load(fullfile(fn.processed,'StandardPd_Spikes'));
InfoTable = readtable(fullfile(fn.processed,sprintf('StandardPd_StimInfo_%s',subject)));
UnitTable = readtable(fullfile(fn.processed,sprintf('UnitDataSummary_%s',subject)));


%%  Prepare table to track units to inspect

if strcmpi(plottype,'comparison')  %could be adapted to others
    
    Temp=table;
    
    SelectUnits = UnitTable(1,1:13);
    SelectUnits.(var1) = nan;
    
    SelectThresh = 80;
    
end

%%

% Set up figures
Figs = struct;
close all

PlotColors = [0.9 0.4 0.4; 0 0 0];
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );
set(0,'DefaultTextInterpreter','none')


scrsz = get(0,'ScreenSize');
figsize = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

% Set axis limits depending on plot type and metric
switch var1
    case 'FR'
        ylimits = [0 100];
    case 'FRnorm'
        ylimits = [-2 15];
    case 'FF'
        ylimits = [0 5];
    case 'PeakTrough'
        ylimits = [0 200];
    case 'PeakTime'
        ylimits = [0 250];
    case 'VS'
        ylimits = [0 0.7];
    case 'FRidx'
        ylimits = [-2 2];
end

if strcmpi(plottype,'comparison')
    xlimits = ylimits;
elseif strcmpi(plottype,'regression')
    switch var2
        case 'FR100ms_idx'
            xlimits = [-2 2];
        case 'FR250ms_idx'
            xlimits = [-2 2];
        case 'prevRate'
            xlimits = [0.5 7];
    end
elseif strcmpi(plottype(1:4),'psth')
    xlimits = [1 250];
    ylimits = [];
end


%%
PlotTypes = {'psth' 'psth_normBase' 'psth_diff' 'comparison' 'regression' 'distribution'};
Var1s = {'FR' 'FRnorm' 'FRidx' 'PeakTrough' 'PeakTime' 'VS' 'VSidx' 'FF'};
Var2s = {'FR100ms_idx' 'FR250ms_idx' 'prevRate'};
Categorizations = {'VS' 'Corr0'};


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
                    
                    % Get data for periodic stimulus
                    T_pdc = filtInfo2(strcmp(filtInfo2.jitter,'0'),:);
                    
                    if isempty(T_pdc), continue, end
                    % Check that just one stimulus
                    if ~numel(unique(T_pdc.trialN))==size(T_pdc,1)
                        keyboard
                    end
                    
                    
                    %% Separate datapoints according to some criterion if desired
                    
                    % Get unit data to filter by response quality
                    UnitData = UnitTable(ismember(UnitTable(:,1:12),T_pdc(1,[1:11 14]),'rows'),:);
                    if size(UnitData,1)>1, keyboard, end
                    
                    dpcat=[];
                    % Designate category
                    switch categorization
                        case 'VS'
                            if UnitData.VS < 0.2
                                dpcat = 'VSbelow'; %below
                            else
                                dpcat = 'VSabove'; %above
                            end
                            
                        case 'Corr'
                            if UnitData.Corr0 < 0
                                dpcat = 'CorrNEG'; %negative correlation
                            elseif isnan(UnitData.Corr0)
                                dpcat = 'Corrins'; %insignificant correlation
                            elseif UnitData.Corr0 >0
                                dpcat = 'CorrPOS'; %negative correlation
                            end
                            
                        case {'all' 'ALL'}
                            dpcat = 'ALL';
                    end
                    if isempty(dpcat), keyboard, end
                    
                    
                    %% Calculate response for periodic stimulus
                    
                    resp_pdc = AmalgamatedSpikes(ismember(InfoTable,T_pdc,'rows'),:);
                    spbin_pdc = binspikecounts(resp_pdc,binsize);
                    FRbin_pdc = mean(spbin_pdc/binsize*1000,1);
                    [FRsmooth_pdc,peakFR_pdc,peakTime_pdc,minFR_pdc] = smoothPSTH(mean(resp_pdc,1));
                    
                    % Skip data points that likely come from junk clusters
                    if mean(FRbin_pdc)<5 || size(resp_pdc,1)<trThresh
                        continue
                    end
                    
                    
                    
                    %% Save FR data for periodic stimulus
                    try
                    if strcmpi(plottype(1:4),'psth')
                        
                        pdcname = [dpcat '_0'];
                        % Create empty matrices for population data if they
                        % don't exist already
                        if ~exist(sprintf('Spikes_%s',pdcname),'var')
                            eval(sprintf('Spikes_%s = [];',pdcname))
                        end
                        % And this datapoint PSTH to big matrix
                        eval(sprintf('Spikes_%s = [Spikes_%s; FRsmooth_pdc];',pdcname,pdcname))
                        
                    end
                    catch
                        aaa=354;
                    end

                    
                    
                    
                    
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
% .........................................................................
% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
                    %% Get unique jitters in this stimulus set
                    
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
                        
                        
                        % Calculate response measures for jitter stimulus
                        
                        resp_jit = AmalgamatedSpikes(ismember(InfoTable,T_jit,'rows'),:);
                        
                        [FRsmooth_jit,peakFR_jit,peakTime_jit,minFR_jit] = smoothPSTH(mean(resp_jit,1));
                        
                        spbin_jit = binspikecounts(resp_jit,binsize);
                        FRbin_jit = mean(spbin_jit/binsize*1000,1);
                        
                        % And calculate jit/pdc indices
                        jitFRidx = (mean(FRbin_jit)-mean(FRbin_pdc)) / mean(FRbin_pdc);%(mean(FRbin_jit)+mean(FRbin_pdc));
                        prev100_idx = (T_jit.prev100msFR(1)-T_pdc.prev100msFR(1)) / T_pdc.prev100msFR(1);%(T_jit.prev100msFR(1)+T_pdc.prev100msFR(1));
                        prev250_idx = (T_jit.prev250msFR(1)-T_pdc.prev250msFR(1)) / T_pdc.prev250msFR(1);%(T_jit.prev250msFR(1)+T_pdc.prev250msFR(1));
                        
                        
                        
                        % Skip data points that likely come from junk clusters
                        if mean(FRbin_jit)<5 || size(resp_jit,1)<trThresh
                            continue
                        end
                        
                        
                        %%   Set up figure handles
                                                
                        %%%%  ALL JITTERS TOGETEHR
                        
                        if ~isfield(Figs,dpcat)
                            
                            Figs.(dpcat) = dpcat;
                            eval(sprintf('%s = figure;',Figs.(dpcat)))
                            set(gcf,'Position',figsize);
                            hold on
                            
                            switch plottype(1:4)
                                case 'comp'
                                    title(sprintf('%s %s during standard period\n%s units',var1,plottype,dpcat))
                                    xlabel(sprintf('%s periodic',var1))
                                    ylabel(sprintf('%s jitter',var1))
                                    plot(xlimits,ylimits,'k','LineWidth',0.5)
                                    set(gca,'xlim',xlimits,'ylim',ylimits)
                                case 'regr'
                                    title(sprintf('%s as a function of %s\n%s units',var1,var2,dpcat))
                                    xlabel(sprintf('%s',var2))
                                    ylabel(sprintf('%s',var1))
                                    plot(xlimits,[0 0],'k','LineWidth',0.5)
                                    plot([0 0],ylimits,'k','LineWidth',0.5)
                                    set(gca,'xlim',xlimits,'ylim',ylimits)
                                case 'psth'
                                    xlim(xlimits)
                                    title(sprintf('%s - during standard 4 Hz period\n%s units',plottype,dpcat))
                                    xlabel('time (ms) in standard pd')
                                    ylabel('Spikes/sec')
                            end
                        end
                        
                        
                        %%%%  EACH JITTER SEPARATLEY
                        
                        ijf = ij{:};
                        ijf(strfind(ijf,'-'))=[];
                        
                        jitfigname = [dpcat '_' ijf];
                        
                        % dont need separate jitter figures for population psth
                        if ~strcmp(plottype,{'psth' 'psth_normBase'})
                            
                        if ~isfield(Figs,jitfigname)
                            
                            Figs.(jitfigname) = jitfigname;
                            eval(sprintf('%s = figure;',Figs.(jitfigname)))
                            set(gcf,'Position',figsize);
                            hold on
                            
                            switch plottype(1:4)
                                case 'comp'
                                    title(sprintf('%s\n%s %s during standard period\n%s units',ij{:},var1,plottype,dpcat))
                                    xlabel(sprintf('%s periodic',var1))
                                    ylabel(sprintf('%s jitter',var1))
                                    plot(xlimits,ylimits,'k','LineWidth',0.5)
                                    set(gca,'xlim',xlimits,'ylim',ylimits)
                                case 'regr'
                                    title(sprintf('%s\n%s as a function of %s\n%s units',ij{:},var1,var2,dpcat))
                                    xlabel(sprintf('%s',var2))
                                    ylabel(sprintf('%s',var1))
                                    plot(xlimits,[0 0],'k','LineWidth',0.5)
                                    plot([0 0],ylimits,'k','LineWidth',0.5)
                                    set(gca,'xlim',xlimits,'ylim',ylimits)
                                case 'psth'
                                    xlim(xlimits)
                                    title(sprintf('%s\n%s - during standard 4 Hz period\n%s units',ij{:},plottype,dpcat))
                                    xlabel('time (ms) in standard pd')
                                    ylabel('Spikes/sec')
                                    
                            end
                            
                        end
                        
                        end %skip for regular population psth
                        
                        
                        
                        %%
                        % Currently comparing distribution of spike counts 
%                         [H,P] = ttest2(spbin_pdc,spbin_jit);
%                         
%                         if sum(H)>4 || any(diff(find(H),2,2)==0) || mean(P,'omitnan')<0.05
%                             sigdp = 1;
%                         else 
%                             sigdp = 2;
%                         end

                                                
%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
                        %% Plot things
%  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
                        try
                            
                            switch plottype(1:4)
                                % * * * * * * * * * * * * * 
                                %    Direct Comparison
                                % * * * * * * * * * * * * * 
                                case 'comp'
                                    
                                    % Get x and y data
                                    switch var1 
                                        case 'FR'
                                            plot_x = mean(FRbin_pdc);
                                            plot_y = mean(FRbin_jit);
                                        case 'FRnorm'
                                            plot_x = (mean(FRbin_pdc) - UnitData.baselineFR)/UnitData.baselineFR;
                                            plot_y = (mean(FRbin_jit) - UnitData.baselineFR)/UnitData.baselineFR;
                                        case 'FF'
                                            keyboard
                                        case 'PeakTime'
                                            plot_x = peakTime_pdc;
                                            plot_y = peakTime_jit;
                                        case 'PeakTrough'
                                            plot_x = peakFR_pdc - minFR_pdc;
                                            plot_y = peakFR_jit - minFR_jit;
                                    end
                                    
                                    % Plot datapoint in everything together figure
                                    eval(sprintf('set(0,''currentfigure'',%s)',Figs.(dpcat)))
                                    hold on
                                    plot(plot_x,plot_y,...
                                        'o','MarkerSize',8,'LineWidth',2,...
                                        'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                                    
                                    % Also plot in corresponding jitter figure
                                    eval(sprintf('set(0,''currentfigure'',%s)',Figs.(jitfigname)))
                                    hold on
                                    plot(plot_x,plot_y,...
                                        'o','MarkerSize',8,'LineWidth',2,...
                                        'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                                    
                                    
                                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    % Pull out units that have a strong response
                                    if strcmpi(var1,'PeakTrough') && plot_y>SelectThresh
                                        Temp.(var1) = plot_y;
                                        SelectUnits = [SelectUnits; [UnitTable(ismember(UnitTable(:,1:12),T_jit(1,[1:11 14]),'rows'),1:13) Temp(1,1)] ];
                                    end
                                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                    
                                    
                                % * * * * * * * * * * * * * 
                                %       Regression
                                % * * * * * * * * * * * * * 
                                case 'regr'
                                    switch var2 
                                        case 'FR100ms_idx'
                                            plot_x = prev100_idx;
                                        case 'FR250ms_idx'
                                            plot_x = prev250_idx;
                                        case 'prevRate'
                                            rv = jitter_LUT(T_jit.CenterRate(1),ij{:});
                                            plot_x = rv(3);
                                            
                                    end
                                    switch var1 %%%%% MORE HERE
                                        case 'FRidx'
                                            plot_y = jitFRidx;
                                    end
                                    
                                    % Plot datapoint in everything together figure
                                    eval(sprintf('set(0,''currentfigure'',%s)',Figs.(dpcat)))
                                    hold on
                                    plot(plot_x,plot_y,...
                                        'o','MarkerSize',8,'LineWidth',2,...
                                        'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                                    
                                    % Also plot in corresponding jitter figure
                                    eval(sprintf('set(0,''currentfigure'',%s)',Figs.(jitfigname)))
                                    hold on
                                    plot(plot_x,plot_y,...
                                        'o','MarkerSize',8,'LineWidth',2,...
                                        'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                                    
                                    
                                % * * * * * * * * * * * * * 
                                %          PSTH
                                % * * * * * * * * * * * * * 
                                case 'psth'
                                    % Create empty matrices for population data if they
                                    % don't exist already
                                    if ~exist(sprintf('Spikes_%s',jitfigname),'var')
                                        eval(sprintf('Spikes_%s = [];',jitfigname))
                                    end
                                    
                                    if strcmpi(plottype,'psth')
                                        % Add this datapoint PSTH to big matrix
                                        eval(sprintf('Spikes_%s = [Spikes_%s; FRsmooth_jit];',jitfigname,jitfigname))
                                        
                                    elseif strcmpi(plottype,'psth_diff')
                                        eval(sprintf('Spikes_%s = [Spikes_%s; FRsmooth_jit-FRsmooth_pdc];',jitfigname,jitfigname))
                                        
                                    end
                            end
                            
                        catch
                            keyboard
                        end
                                                
                        
                    end
                    
                    %%
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Pull out units that have a strong response
                    if strcmpi(var1,'PeakTrough') && plot_x>SelectThresh  && strcmpi(plottype(1:4),'comp')
                        Temp.(var1) = plot_y;
                        SelectUnits = [SelectUnits; [UnitData(1,1:13) Temp(1,1) ]];
                    end
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                end
            end
        end
    end
end


%% Plot population PSTHs

if  strcmpi(plottype(1:4),'psth')
    
    C = who('-regexp', ['Spikes_' categorization]);
    
    cats = unique(strtok(extractAfter(C,'Spikes_'),'_'));
    
    for icat = cats'
        
        C = who('-regexp', ['Spikes_' icat{:}]);
        
        
        % Go through each jitter
        for icj = C' 
            
            % Require at least 5 units for this jitter
            if size(eval(icj{:}),1)<5 && ~onlySU
                continue
            end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %  Regular population PSTH
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if strcmpi(plottype,'psth')
                
                eval(sprintf('figure(%s);',icat{:}))
                hold on
                
                if strcmp( strtok(extractAfter(icj{:},[icat{:} '_']),'_'), '0' )
                    
                    fill([1:250 250:-1:1],...
                        [mean(eval(icj{:}),1)+std(eval(icj{:}),1) fliplr(mean(eval(icj{:}),1)-std(eval(icj{:}),1))],...
                        [0 0 0],'FaceAlpha',0.4)
                    
                else
                    plot( mean(eval(icj{:}),1)' ,...
                        'LineWidth',2,...
                        'Color', ALLcolors( strcmp( strtok(extractAfter(icj{:},[icat{:} '_']),'_'), strtok(plotOptions.colSelect,'_') ), : ));
                end
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %  Jitter - Periodic PSTHs
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            elseif strcmpi(plottype,'psth_diff')
                
                if strcmp( strtok(extractAfter(icj{:},[icat{:} '_']),'_'), '0' )
                    continue
                end
                
                % Plot each difference in the appropriate jitter figure
                eval( sprintf('figure(%s_%s);',icat{:},extractAfter(icj{:},[icat{:} '_'])) )
                hold on
                
                plot(repmat(1:250,[size(eval(icj{:}),1) 1])' , eval(icj{:})' ,...
                    'Color', ALLcolors( strcmp( strtok(extractAfter(icj{:},[icat{:} '_']),'_'), strtok(plotOptions.colSelect,'_') ), : ),...
                    'LineWidth',2);
                
                % Plot mean difference in plot all together
                eval(sprintf('figure(%s);',icat{:}))
                hold on
                
                plot( mean(eval(icj{:}),1)' ,...
                    'LineWidth',2,...
                    'Color', ALLcolors( strcmp( strtok(extractAfter(icj{:},[icat{:} '_']),'_'), strtok(plotOptions.colSelect,'_') ), : ));
                
            end
        end
        
    end
    
end


%% Save figures


% Set save directory
if onlySU
    savedir = fullfile(fn.standardPd,'onlySU');
else
    savedir = fn.standardPd;
end

savedir = fullfile(savedir,plottype);

if ~isempty(var1)
    savedir = fullfile(savedir,var1);
end
if ~isempty(var2)
    savedir = fullfile(savedir,var2);
end

if ~exist(savedir,'dir')
    mkdir(savedir)
end


% Save figures
FigNames = fieldnames(Figs);

for ifn = FigNames'
    
    [categ,jit] = strtok(ifn{1},'_');
    
    if isempty(jit)
        jit = '_AllJits';
    end
    
    switch plottype(1:4)
        case 'comp'
            savename   = [var1 '_' plottype '_' categ jit];
        case 'regr'
            savename   = [var1 '_' var2 '_' plottype '_' categ jit];
        case 'psth'
            savename   = [plottype '_' categ jit];
    end
    
    eval(sprintf('figure(%s)',Figs.(ifn{1})))
    
    eval(sprintf('set(%s,''PaperOrientation'',''landscape'')', Figs.(ifn{1}) ))
    print(eval(Figs.(ifn{1})),'-dpdf',fullfile(savedir,savename),'-bestfit');
    
end



%% Save SelectUnits table (if applicable)

if exist('SelectUnits','var')
    
    if isnan(SelectUnits(1,:).(var1))
        SelectUnits(1,:) = [];
    end
    
    tablesavename = sprintf('SelectUnits_%s_%i',var1,SelectThresh);
    writetable(SelectUnits,fullfile(savedir,tablesavename));
    
end




end





