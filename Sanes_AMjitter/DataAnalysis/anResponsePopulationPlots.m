function anResponsePopulationPlots(subject,METRICS)



clear DataTable
global DataTable
% global session unit

fn = set_paths_directories;


%%
if nargin<2 && ~exist('METRICS','var')
%****************************************%
METRICS = {'FR' 'FF' 'FFavPds' 'VS' 'RS' 'standardFR' 'maxCorr' 'shftCorr'};
%****************************************%
end

%% Table info

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
DataTable.baselineFR   = 0;
DataTable.FR           = 0;
DataTable.FF           = 0;
DataTable.FF_avPds     = 0;
DataTable.VS           = 0;
DataTable.RS           = 0;
DataTable.standardFR   = 0;
DataTable.Corr_max     = 0;
DataTable.Corr_shft    = 0;
DataTable.DiscrmVar    = 'placeholder';


%% FOR EACH INDEPENDENT VARIABLE

for iV = { 'depth' 'jitter' }

    indVar = iV{1};
    nmfieldname = ['nm_' indVar];
    close all


%% Plot info
Figs = struct;
[~, plotOptions] = setOptions;
plotOptions.plotThresh = false;
switch indVar
    case 'depth'
        plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
        ALLcolors = copper( numel(plotOptions.colSelect) );
        xlimits = [-30 0];
    case 'jitter'
%         ALLcolors = winter( numel(plotOptions.colSelect) );
        ALLcolors = [ 0   1 50  ;...
                      0  11 90  ;...
                      0  43 45  ;...
                      1  75  0 ] ./ 100;
        plotOptions.xLabel = 'jittered from periodic AM (x/100 = range in log2 units around middle rate)';
        xlimits = [0 250];
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
        
        % *****  SKIP IF << create_Data_struct >> (new version) NOT RUN YET  ***** %
        if ~isfield(Data.(unit{:}),'stimdata'), continue, end
        % ************************************************************************ %
        
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
        
        % Otherwise, plot that data
        for is = get_sess_data
            for ip = 1:numel(stimdata(is).pars)                
                
                % Get stimulus parameters and behavioral state
                stimpars = stimdata(is).pars(ip).pars;
                behavs = unique(stimdata(is).pars(ip).stimvals(:,3));
                
                for ib = 1:numel(behavs)  %separate behavioral states
                    
                    for im = 1:numel(METRICS)
                        
                        this_metric = METRICS{im};
                        
                        %% Set up figure handles
                        
                        parname = sprintf('%s_%i_%i_%i_%i_%s',this_metric,stimpars(1),stimpars(2),stimpars(3),stimpars(4),behavs(ib));
                        
                        if ~isfield(Figs,parname)
                            
                            % * * * * * * * * * * * * * * * * * * * * * * * *
                            % Analysis Figures
                            % * * * * * * * * * * * * * * * * * * * * * * * *
                            Figs.(parname).anF = ['anF' num2str(numel(fieldnames(Figs)))];
                            eval(sprintf('%s = figure;',Figs.(parname).anF))
                            
                            scrsz = get(0,'ScreenSize');
                            eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).anF))
                            
                            
                            % * * * * * * * * * * * * * * * * * * * * * * * *
                            % Jitter Condition Difference Figures
                            % * * * * * * * * * * * * * * * * * * * * * * * *
                            if strcmp(indVar,'depth')
                                Figs.(parname).DanF = ['DanF' num2str(numel(fieldnames(Figs)))];
                                eval(sprintf('%s = figure;',Figs.(parname).DanF))
                                
                                scrsz = get(0,'ScreenSize');
                                eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).DanF))
                            end
                            
                            % * * * * * * * * * * * * * * * * * * * * * * * *
                            % Histogram Figure
                            % * * * * * * * * * * * * * * * * * * * * * * * *
%                             Figs.(parname).histF = ['histF' num2str(numel(fieldnames(Figs)))];
%                             eval(sprintf('%s = figure;',Figs.(parname).histF))
%                             eval(sprintf('set(%s,''Position'',[1 scrsz(4) scrsz(3) scrsz(4)],''Nextplot'',''add'');',Figs.(parname).histF))
                            
                        end
                        
                        
                        %% Get response analysis data
                        
                        rasters = stimdata(is).pars(ip).(nmfieldname);
                        rasters = rasters(strcmp({rasters.behaving},behavs(ib)));
                        
                        % Skip the incomplete depth functions (due to not
                        % enough trials while drinking
                        if     numel(rasters)<9 && strcmp(behavs(ib),'D') && strcmp(indVar,'depth')
                            continue
                        elseif numel(rasters)<9 && strcmp(behavs(ib),'P') && strcmp(indVar,'depth')
                            keyboard
                        end
                        
                        % Plot
                        eval(sprintf('figure(%s)',Figs.(parname).anF))
                        hold on
                        dp = an_responses(this_metric,rasters,indVar,subject,plotOptions,ALLcolors,xlimits,0);
                        % * * * * next: add response data to summary table
                        
                        xlabel(plotOptions.xLabel)
                        titlestr = sprintf('%s response functions\n%s', this_metric, parname);
                        title(titlestr)
                        
                        pause(0.1)
                        
                        % Plot jitter condition difference
                        if strcmp(indVar,'depth')
                            
                            eval(sprintf('figure(%s)',Figs.(parname).DanF))
                            hold on
                            an_responses(this_metric,rasters,indVar,subject,plotOptions,ALLcolors,xlimits,1);
                            
                            xlabel(plotOptions.xLabel)
                            titlestr = sprintf('%s responses, diff between conditions\n%s', this_metric, parname);
                            title(titlestr)
                            
                        end
                        
                        pause(0.1)
                        
%                         %% Add row to Data Table
%                         
%                         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         add_datapoints_table(stimdata(is).pars(ip).stimvals,...
%                             clustruct(is).pars(ip).(nmfieldname),...
%                             session, unit{:}, Data.(unit{:}).labels(1,2),...
%                             stimpars, behavs(ib), indVar )
%                         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                         
                        
                    end % METRICS
                end %behavs
            end %pars
        end %stim block
    end %unit
    
end %session

% Remove first row of Data Table if it's still the placeholder
% if strcmp(DataTable{1,1},'placeholder')
%     DataTable(1,:)=[];
% end

%% Histogram plots

% % Filter the table data according to basic parameters
% DT = DataTable(strcmp(DataTable.DiscrmVar,indVar),:);
% Pars = unique( [DT.Noiseband DT.dBSPL DT.CenterRate ] ,'rows');
% Behs = unique(DT.BehavState);
% 
% for ip = 1:size(Pars,1)
%     for ib = Behs'
%         
%         FilteredTable = DT( DT.BehavState==ib & ...
%                                all(DT.Noiseband == [Pars(ip,1) Pars(ip,2)],2) & ...
%                                    DT.dBSPL == Pars(ip,3) & ...
%                                    DT.CenterRate == Pars(ip,4) , : );
%         
%         % Get independent and condition variable values
%         indVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,indVar)},'stable');
%         indVar_vals(strcmp(indVar_vals,'0')) = [];
%         switch indVar
%             case 'jitter'
%                 condVar = 'depth';
%                 condVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,condVar)});
%                 condVar_vals(strcmp(condVar_vals,'0')) = [];
%             case 'depth'
%                 condVar = 'jitter';
%                 condVar_vals = unique(FilteredTable{:,strcmpi(DT.Properties.VariableNames,condVar)});
%         end
%         
%         
%         % Prepare figure
%         beh=char(ib);
%         parname = sprintf('pars_%i_%i_%i_%i_%s',Pars(ip,1),Pars(ip,2),Pars(ip,3),Pars(ip,4),beh(1));
%         
%         if ~isfield(Figs,parname), continue, end
%         
%         eval(sprintf('figure(%s)',Figs.(parname).histF))
%         clf
%         hold on
%         ih=0;
%         for isp = 1:numel(indVar_vals)
%             
%             % Plot this subplot
%             sp(isp) = subplot(numel(indVar_vals),1,isp);
%             hold on
%             for ic = 1:numel(condVar_vals)
%                 
%                 ih = ih+1; 
%                 ftd = (strcmp( FilteredTable{:, strcmpi(DT.Properties.VariableNames,indVar)},  indVar_vals(isp) ) )...
%                     & (strcmp( FilteredTable{:, strcmpi(DT.Properties.VariableNames,condVar)}, condVar_vals(ic) ) );
%                 
%                 hp(ih) = histogram(FilteredTable.(selectCol)(ftd));
%                 hp(ih).Normalization = 'count';
%                 hp(ih).BinWidth = 0.1;
%                 hp(ih).FaceColor = ALLcolors( strcmp(strtok(condVar_vals(ic),'_'),strtok(plotOptions.colSelect,'_')) ,:);
%                 hp(ih).EdgeColor = ALLcolors( strcmp(strtok(condVar_vals(ic),'_'),strtok(plotOptions.colSelect,'_')) ,:);
%                                 
%             end
%             
%             title(indVar_vals(isp),'FontName',plotOptions.fontName,'FontSize', plotOptions.labelSize)
%             if isp==numel(indVar_vals)
%                 xlabel(['d'' values from ' which_classifier],'FontName',plotOptions.fontName,'FontSize', plotOptions.labelSize)
%             end
%         end
%         
%         % Finish figure
%         set(gca,'xlim',[0 ceil(2*max(cell2mat(cellfun(@max,{hp.Data},'UniformOutput',false))))/2])
%         linkaxes(sp,'xy')
%         
%         titlestr = sprintf('Histogram of all units'' dprime values, discriminating %s\n%s\n%s', indVar, parname, which_classifier);
%         suptitle(titlestr)
%         
%         clear hp sp
%     end
% end

%% Finish and save figures

savedir = fullfile(fn.processed,'^Population','ResponseMetrics',indVar);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

for ifn = fieldnames(Figs)'
    
    parcell = strsplit(ifn{1},'_');
    
    for ift = fieldnames(Figs.(ifn{1}))'
        
%         % Set y axis limits
%         eval(sprintf('figure(%s)',Figs.(ifn{1}).anF))
%         ylms = get(gca,'ylim');
        
        if strcmp(ift{1}(1),'D')
            savename   = sprintf('AM%sHz_%s-%s_%sdB_%s_%s_allUnits_%s-CondDiff',...
                parcell{5},parcell{2},parcell{3},parcell{4},parcell{6},subject,parcell{1});
        else
            savename   = sprintf('AM%sHz_%s-%s_%sdB_%s_%s_allUnits_%s',...
                parcell{5},parcell{2},parcell{3},parcell{4},parcell{6},subject,parcell{1});
        end
        
        % Save figure
        set(gcf,'PaperOrientation','landscape');
        eval(sprintf('print(%s,''-dpdf'',fullfile(savedir,savename),''-bestfit'');',Figs.(ifn{1}).(ift{1})))
        print(eval(Figs.(ifn{1}).(ift{1})),'-dpdf',fullfile(savedir,savename),'-bestfit');
        
    end
    
end


end % for indVar



%% Save the data table
% tablefn = sprintf('UnitDataSummary_%s',subject);
% writetable(DataTable,fullfile(fn.processed,tablefn));


end









