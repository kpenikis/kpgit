function updatePopulationPlots(subject,which_classifier)


close all
clear DataTable
global DataTable

%% Set up useful variables

%****************************************%
%****************************************%
% IMPORTANT DATA CHOICES

NM_DIFF = 1;

if nargin<2
%     which_classifier = 'classifFR0';
%     which_classifier = 'classifSpV100';
%     which_classifier = 'classifSpV25';
    which_classifier = 'classifSpV10';
%     which_classifier = 'classifCorr0';
%     which_classifier = 'formulaFR';  % must re-run getALL_nmdata and then this will work (dont forget to change data_table function too
end
%****************************************%
%****************************************%


%% Table info
fn = set_paths_directories;

DataTable = table;
%'Session' 'cluname' 'SU/MU' 'HP LP' 'dBSPL' 'CenterRate' 'BehavState' 'Depth' 'Jitter' 'dprime' 'indVar'
DataTable.Session      = 'placeholder';
DataTable.cluname      = 'placeholder';
DataTable.UnitType     = discretize(4,[0.5 1.5 2.5 3.5 4.5],'categorical',{'unk', 'SU', 'MU', 'noise'});
DataTable.Noiseband    = [0 0];
DataTable.dBSPL        = 0;
DataTable.CenterRate   = 0;
DataTable.BehavState   = categorical({'P'},{'P' 'D' 'A'},{'Passive' 'Drinking' 'Active Behavior'});
DataTable.depth        = 0;
DataTable.jitter       = 'placeholder';
DataTable.dprimeFR_cls = 0;
DataTable.dprime100ms  = 0;
DataTable.dprime25ms   = 0;
DataTable.dprime10ms   = 0;
DataTable.dprimeCorr   = 0;
DataTable.dprimeFR_fml = 0;
DataTable.DiscrmVar    = 'placeholder';

% Set Variable Name corresponding to desired dprime method
switch which_classifier
    case 'classifFR0'
        selectCol = 'dprimeFR_cls';
        METRIC = 'FR';
    case 'classifSpV100'
        selectCol = 'dprime100ms';
        METRIC = 'FR';
    case 'classifSpV25'
        selectCol = 'dprime25ms';
        METRIC = 'FR';
    case 'classifSpV10'
        selectCol = 'dprime10ms';
        METRIC = 'FR';
    case 'classifCorr0'
        selectCol = 'dprimeCorr';
        METRIC = 'Corr';
end
    

%% FOR EACH INDEPENDENT VARIABLE

for iV = {'depth' 'jitter'}

    indVar = iV{1};
    nmfieldname = ['nm_' indVar];


%% Plot info
Figs = struct;
[~, plotOptions] = setOptions;
plotOptions.plotThresh = false;
switch indVar
    case 'depth'
        plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
        ALLcolors = copper( numel(plotOptions.colSelect) );
        xlimits = [-20 0];
        ydplim  = [0 3];
    case 'jitter'
%         ALLcolors = winter( numel(plotOptions.colSelect) );
        ALLcolors = [ 0   1 50  ;...
                      8  26 40  ;...
                     15  51 31  ;...
                     23  75 22  ;...
                     31 100 12 ] ./ 100;
        plotOptions.xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
        xlimits = [0 250];
        ydplim  = [0 2];
end
set(0,'DefaultTextInterpreter','none')
set(0,'defaultaxesfontsize',plotOptions.labelSize)


%% Get all processed Data files
searchname = sprintf('%s_sess-*_Data.mat',subject);
alldatafiles = dir(fullfile(fn.processed,subject,searchname));

for ii = 1:numel(alldatafiles)
    
    % Load Data file
    clear Data
    load(fullfile(fn.processed,subject,alldatafiles(ii).name));
    
    % Get this session label and set directories
    [~,session] = strtok(alldatafiles(ii).name,'-');
    session = session(2:3);
    fn = set_paths_directories(subject,session);
    
    % Load each cluster file
    allclusters = fieldnames(Data);
    allclusters = allclusters(strncmp('ch',allclusters,2));
    
    
    for unit = allclusters'
        
        % *****  SKIP IF << create_Data_struct >> NOT RUN YET  ***** %
        if ~isfield(Data.(unit{:}),'stimdata'), continue, end
        % ********************************************************** %
        
        % Check if this session contains the discriminated variable
        stimdata = Data.(unit{:}).stimdata;
        get_sess_data = [];
        for is = 1:numel(stimdata)
            if isfield(stimdata(is).pars,nmfieldname)
                get_sess_data = [get_sess_data is];
            end
        end
        
        if isempty(get_sess_data), continue, end
        
        %% If there is relevant data in this session, load cluster data file
        
        clear clustruct
        cluefilename = sprintf( '%s_sess-%s_%s',subject,session,unit{:});
        clustruct = load(fullfile(fn.sess_data,cluefilename));
        clustruct = clustruct.(unit{:});
        
        
        % *****  SKIP IF << getALL_nmData >> NOT RUN YET  ***** %
        if ~isfield(clustruct,'pars'), continue, end
        % ***************************************************** %
        
        % Otherwise, plot that data
        for is = get_sess_data
            for ip = 1:numel(clustruct(is).pars)                
                
                % Get stimulus parameters and behavioral state
                stimpars = stimdata(is).pars(ip).pars;
                behavs = unique(stimdata(is).pars(ip).stimvals(:,3));
                
                for ib = 1:numel(behavs)  %separate behavioral states
                    
                    %% Set up figure handles
                    
                    parname = sprintf('pars_%i_%i_%i_%i_%s',stimpars(1),stimpars(2),stimpars(3),stimpars(4),behavs(ib));
                    
                    if ~isfield(Figs,parname)
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        % Neurometric Figures
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        Figs.(parname).nmF = ['nmF' num2str(numel(fieldnames(Figs))+1)];
                        eval(sprintf('%s = figure;',Figs.(parname).nmF))
                        
                        scrsz = get(0,'ScreenSize');
                        eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).nmF))
                        
                        nsubplots = 1;
                        for isp = 1:nsubplots
                            Figs.(parname).nmS = ['nmS' num2str(numel(fieldnames(Figs)))];
                            eval(sprintf('%s=subplot(1,%i,%i);',Figs.(parname).nmS,nsubplots,isp))
                        end
                        
                        % Neurometric Difference Figures
                        if NM_DIFF && strcmp(indVar,'depth')
                            Figs.(parname).nmdF = ['nmdF' num2str(numel(fieldnames(Figs)))];
                            eval(sprintf('%s = figure;',Figs.(parname).nmdF))
                            
                            scrsz = get(0,'ScreenSize');
                            eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).nmdF))
                        end
                        
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        % Analysis Figure
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        Figs.(parname).anF = ['anF' num2str(numel(fieldnames(Figs)))];
                        eval(sprintf('%s = figure;',Figs.(parname).anF))
                        
                        scrsz = get(0,'ScreenSize');
                        eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).anF))
                        
                        nsubplots = 1;
                        for isp = 1:nsubplots
                            Figs.(parname).anS = ['nmS' num2str(numel(fieldnames(Figs)))];
                            eval(sprintf('%s=subplot(1,%i,%i);',Figs.(parname).nmS,nsubplots,isp))
                        end
                        
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        % Histogram Figure
                        % * * * * * * * * * * * * * * * * * * * * * * * *
                        Figs.(parname).histF = ['histF' num2str(numel(fieldnames(Figs)))];
                        eval(sprintf('%s = figure;',Figs.(parname).histF))
                        eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).histF))
                        
                    end
                    
                    
                    %% Plot neurometric data
                    
                    % Get data to plot
                    this_behav = strcmp(behavs(ib),{clustruct(is).pars(ip).(nmfieldname).(which_classifier).behav});
                    dp_struct = clustruct(is).pars(ip).(nmfieldname).(which_classifier)(this_behav).dprime;
                    
                    % Plot the data
                    if NM_DIFF && strcmp(indVar,'depth')
                        
                        % Set current figure & subplot
                        eval(sprintf('figure(%s)',Figs.(parname).nmdF))
                        hold on
                        
                        population_plot_nmdifference(dp_struct,...
                            clustruct(is).pars(ip).(nmfieldname).stim,...
                            plotOptions, ALLcolors, [], xlimits)
                        
                        set(gca,'ylim',[-0.75 0.75])
                        titlestr = sprintf('Difference between conditions: neurometric %s functions\n%s\n%s', indVar, parname, which_classifier);
                        suptitle(titlestr)
                        
                    end
                    
                    % Set current figure & subplot
                    eval(sprintf('figure(%s)',Figs.(parname).nmF))
                    hold on
                    eval(sprintf('hS = subplot(%s);',Figs.(parname).nmS))
                    
                    population_plot_neurometrics(dp_struct,...
                        clustruct(is).pars(ip).(nmfieldname).stim,...
                        plotOptions, ALLcolors, hS, xlimits)
                    
                    set(gca,'ylim',ydplim)
                    titlestr = sprintf('Neurometric %s functions\n%s\n%s', indVar, parname, which_classifier);                    
                    suptitle(titlestr)
                    
                    
                    %% Get response analysis data
                    
                    rasters = stimdata(is).pars(ip).(nmfieldname);
                    rasters = rasters(strcmp({rasters.behaving},behavs(ib)));
                    
                    eval(sprintf('figure(%s)',Figs.(parname).anF))
                    hold on
                    dp = an_responses(METRIC,rasters,indVar,subject,plotOptions,ALLcolors,xlimits);
                    
                    % * * * * next task: add response data to data table
                    
                    xlabel(plotOptions.xLabel)
                    titlestr = sprintf('%s response functions\n%s', METRIC, parname);
                    title(titlestr)
                    
                    %% Add row to Data Table
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    add_datapoints_table(stimdata(is).pars(ip).stimvals,...
                        clustruct(is).pars(ip).(nmfieldname),...
                        session, unit{:}, Data.(unit{:}).labels(1,2),...
                        stimpars, behavs(ib), indVar )
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    
                end %behavs
            end %pars
        end %stim block
    end %unit
    
end %session

% Remove first row of Data Table if it's still the placeholder
if strcmp(DataTable{1,1},'placeholder')
    DataTable(1,:)=[];
end

%% Histogram plots


% Filter the table data according to basic parameters
DT = DataTable(strcmp(DataTable.DiscrmVar,indVar),:);
Pars = unique( [DT.Noiseband DT.dBSPL DT.CenterRate ] ,'rows');
Behs = unique(DT.BehavState);

for ip = 1:size(Pars,1)
    for ib = Behs'
        
        FilteredTable = DT( DT.BehavState==ib & ...
                               all(DT.Noiseband == [Pars(ip,1) Pars(ip,2)],2) & ...
                                   DT.dBSPL == Pars(ip,3) & ...
                                   DT.CenterRate == Pars(ip,4) , : );
        
        % Get independent and condition variable values
        indVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,indVar)},'stable');
        indVar_vals(strcmp(indVar_vals,'0')) = [];
        switch indVar
            case 'jitter'
                condVar = 'depth';
                condVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,condVar)});
                condVar_vals(strcmp(condVar_vals,'0')) = [];
            case 'depth'
                condVar = 'jitter';
                condVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,condVar)});
        end
        
        
        % Prepare figure
        beh=char(ib);
        parname = sprintf('pars_%i_%i_%i_%i_%s',Pars(ip,1),Pars(ip,2),Pars(ip,3),Pars(ip,4),beh(1));
        
        if ~isfield(Figs,parname), continue, end
        
        eval(sprintf('figure(%s)',Figs.(parname).histF))
        clf
        hold on
        ih=0;
        for isp = 1:numel(indVar_vals)
            
            % Plot this subplot
            sp(isp) = subplot(numel(indVar_vals),1,isp);
            hold on
            for ic = 1:numel(condVar_vals)
                
                ih = ih+1; 
                ftd = (strcmp( FilteredTable{:, strcmpi(DT.Properties.VariableNames,indVar)},  indVar_vals(isp) ) )...
                    & (strcmp( FilteredTable{:, strcmpi(DT.Properties.VariableNames,condVar)}, condVar_vals(ic) ) );
                
                hp(ih) = histogram(FilteredTable.(selectCol)(ftd));
                hp(ih).Normalization = 'count';
                hp(ih).BinWidth = 0.1;
                hp(ih).FaceColor = ALLcolors( strcmp(strtok(condVar_vals(ic),'_'),strtok(plotOptions.colSelect,'_')) ,:);
                hp(ih).EdgeColor = ALLcolors( strcmp(strtok(condVar_vals(ic),'_'),strtok(plotOptions.colSelect,'_')) ,:);
                                
            end
            
            title(indVar_vals(isp),'FontName',plotOptions.fontName,'FontSize', plotOptions.labelSize)
            if isp==numel(indVar_vals)
                xlabel(['d'' values from ' which_classifier],'FontName',plotOptions.fontName,'FontSize', plotOptions.labelSize)
            end
        end
        
        % Finish figure
        set(gca,'xlim',[0 ceil(2*max(cell2mat(cellfun(@max,{hp.Data},'UniformOutput',false))))/2])
        linkaxes(sp,'xy')
        
        titlestr = sprintf('Histogram of all units'' dprime values, discriminating %s\n%s\n%s', indVar, parname, which_classifier);
        suptitle(titlestr)
        
        clear hp sp
    end
end



%% Save figures




end % for indVar

keyboard

%% Save the data table
tablefn = sprintf('UnitDataSummary_%s',subject);
writetable(DataTable,fullfile(fn.processed,tablefn));


end









