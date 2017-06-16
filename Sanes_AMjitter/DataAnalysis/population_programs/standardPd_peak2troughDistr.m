function TotalUnitCount = standardPd_peak2troughDistr(subject,CRITERION)

% Next:
%  - normalize by baseline rate
%  - periodic boxplot in blue
%  - separate drinking and passive

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
set(0,'DefaultAxesFontSize',16)
set(0,'defaulttextfontsize',16)
set(0,'DefaultTextInterpreter','none')

close all
PlotColors = [0.9 0.4 0.4; 0 0 0];
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );

alljitters = unique(UnitTable.jitter);
alljitters = sort_jitters(alljitters);

if strcmpi(CRITERION,'all')
    
    % all datapoints
    hf(1) = figure; hold on
    title('peak FR - trough FR, ALL datapoints')
    ylabel('peak FR - trough FR')
    
    hf_pk(1) = figure; hold on
    title('peak FR, ALL datapoints')
    ylabel('peak FR')
    
elseif strcmpi(CRITERION,{'VS' 'Corr'})
    % below Crit
    hf(1) = figure; hold on
    title(['peak FR - trough FR, below ' CRITERION ' crit'])
    ylabel('peak FR - trough FR')
    
    hf_pk(1) = figure; hold on
    title(['peak FR, below ' CRITERION ' crit'])
    ylabel('peak FR')
    
    % above Crit
    hf(2) = figure; hold on
    title(['peak FR - trough FR, above ' CRITERION ' crit'])
    ylabel('peak FR - trough FR')
    
    hf_pk(2) = figure; hold on
    title(['peak FR, above ' CRITERION ' crit'])
    ylabel('peak FR')
    
elseif strcmpi(CRITERION,'BehState')
    
    % Passive
    hf(1) = figure; hold on
    title(['peak FR - trough FR, below ' CRITERION ' crit'])
    ylabel('peak FR - trough FR')
    
    hf_pk(1) = figure; hold on
    title(['peak FR, below ' CRITERION ' crit'])
    ylabel('peak FR')
    
    % Drinking
    hf(2) = figure; hold on
    title(['peak FR - trough FR, above ' CRITERION ' crit'])
    ylabel('peak FR - trough FR')
    
    hf_pk(2) = figure; hold on
    title(['peak FR, above ' CRITERION ' crit'])
    ylabel('peak FR')
        
end



%% Go through data

TotalUnitCount = 0;

sessions = unique(InfoTable.Session,'stable');
for isess = sessions'
    for iu = unique(InfoTable(strcmp(InfoTable.Session,isess{:}),:).cluname,'stable')'
        
        if onlySU && ~any(strcmp(InfoTable(strcmp(InfoTable.Session,isess{:}) & strcmp(InfoTable.cluname,iu{:}),:).UnitType,'SU'))
            continue
        end
        
        TotalUnitCount = TotalUnitCount+1;
        
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
                        resp_jit  = AmalgamatedSpikes(ismember(InfoTable,T_jit,'rows'),:);
                        spbin_jit = binspikecounts(resp_jit,binsize);
                        FRbin_jit = mean(spbin_jit/binsize*1000,1);
                        
                        % Skip data points that likely come from junk clusters
                        if mean(FRbin_jit)<5 || mean(FRbin_pdc)<5
                            continue
                        end
                        
                        %% Plot things
                        
                        jit = find(strcmp(ij,alljitters));

                        figure(hf(dpcat));
                        plot( jit, max(FRbin_jit)-min(FRbin_jit),...
                            'o','MarkerSize',16,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                        
                        figure(hf_pk(dpcat));
                        plot( jit, max(FRbin_jit),...
                            'o','MarkerSize',16,'LineWidth',2,...
                            'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ));
                        
                        % If this is the first jitter of the datapoint,
                        % also plot the periodic data
                        if find(strcmp(ij,jitters))==1
                            figure(hf(dpcat));
                            plot( 1, max(FRbin_pdc) - min(FRbin_pdc),...
                                'o','MarkerSize',16,'LineWidth',2,...
                                'Color', ALLcolors( strcmp('0',strtok(plotOptions.colSelect,'_')), : ));
                            figure(hf_pk(dpcat));
                            plot( 1, max(FRbin_pdc),...
                                'o','MarkerSize',16,'LineWidth',2,...
                                'Color', ALLcolors( strcmp('0',strtok(plotOptions.colSelect,'_')), : ));
                        end
                        
                    end
                end
            end
        end
    end
end


%% Get data from each plot and make box figures
for idp = 1:numel(hf)
    
    XData=[]; YData=[]; Lines=[];
    
    % Peak - trough
    Lines = get(get(hf(idp),'Children'),'Children');
    XData = [Lines.XData];
    YData = [Lines.YData];
    
    figure(hf(idp))
    bp=boxplot(YData,XData,'PlotStyle','traditional','Jitter',0,...
        'Labels',alljitters,'LabelOrientation','inline','Widths',0.75);
    set(bp,'LineWidth',2,'Color',[0 0 0])
    
    XData=[]; YData=[]; Lines=[];
    
    % Peak only
    Lines = get(get(hf_pk(idp),'Children'),'Children');
    XData = [Lines.XData];
    YData = [Lines.YData];
    
    figure(hf_pk(idp))
    bp=boxplot(YData,XData,'PlotStyle','traditional','Jitter',0,...
        'Labels',alljitters,'LabelOrientation','inline','Widths',0.75);
    set(bp,'LineWidth',2,'Color',[0 0 0])
    
end


%% Save figures

savedir = fullfile(fn.standardPd,'Peak-Trough');
if onlySU
    savedir = fullfile(savedir,'onlySU');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end

if strcmpi(CRITERION,'all')

% hf
set(hf(1),'PaperOrientation','landscape');
print(hf(1),fullfile(savedir,'PeakTroughDiff_ALL'),'-dpdf','-bestfit');

% hf_pk
set(hf_pk(1),'PaperOrientation','landscape');
print(hf_pk(1),fullfile(savedir,'PeakFR_ALL'),'-dpdf','-bestfit');

else     % categorized
    
% hf
set(hf(1),'PaperOrientation','landscape');
print(hf(1),fullfile(savedir,['PeakTroughDiff_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf(2),'PaperOrientation','landscape');
print(hf(2),fullfile(savedir,['PeakTroughDiff_aboveCrit' CRITERION]),'-dpdf','-bestfit');

% hf_pk
set(hf_pk(1),'PaperOrientation','landscape');
print(hf_pk(1),fullfile(savedir,['PeakFR_belowCrit' CRITERION]),'-dpdf','-bestfit');

set(hf_pk(2),'PaperOrientation','landscape');
print(hf_pk(2),fullfile(savedir,['PeakFR_aboveCrit' CRITERION]),'-dpdf','-bestfit');

end



end


