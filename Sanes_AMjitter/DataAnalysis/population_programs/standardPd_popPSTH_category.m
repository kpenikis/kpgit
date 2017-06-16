function standardPd_popPSTH_category(subject)

global fn 

fn = set_paths_directories;

% Load data
load(fullfile(fn.processed,'StandardPd_Spikes'));
InfoTable = readtable(fullfile(fn.processed,sprintf('StandardPd_StimInfo_%s',subject)));
UnitTable = readtable(fullfile(fn.processed,sprintf('UnitDataSummary_%s',subject)));


%%

binsize = 25;

% Set up figures
close all
PlotColors = [0.9 0.4 0.4; 0 0 0];
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );

plot_order = [10  2  6 16 ;...
               9  3  7 15 ;...
              12  1 14  5 ;...
               4 11  8 13 ];
plot_order = flipud(rot90(plot_order));

% Population average - negative Corr
hf_pop_neg = figure; hold on
xlim([1 250/binsize])
% ylim([0 60])
title(sprintf('PSTH of activity during standard 4 Hz period\nfor all units with negative Corr to Sound waveform (no shift)'))
xlabel('time (ms) in standard pd')
ylabel('FR')

% Population average - insign. Corr
hf_pop_ins = figure; hold on
xlim([1 250/binsize])
% ylim([0 80])
title(sprintf('PSTH of activity during standard 4 Hz period\nfor all units with insign. Corr to Sound waveform (no shift)'))
xlabel('time (ms) in standard pd')
ylabel('FR')

% Population average - positive Corr
hf_pop_pos = figure; hold on
xlim([1 250/binsize])
% ylim([0 80])
title(sprintf('PSTH of activity during standard 4 Hz period\nfor all units with positive Corr to Sound waveform (no shift)'))
xlabel('time (ms) in standard pd')
ylabel('FR')


% Go through data
is = 0;
sessions = unique(InfoTable.Session,'stable');
for isess = sessions'
    
    if ~any(strcmp(isess,unique(UnitTable.Session)))
        warning('Unit data has not yet been collected for this session. Re-run <<anResponsePopulationPlots>>')
        keyboard
        continue
    end
    
    is=is+1;
    hfp(is) = figure; hold on
    suptitle(['session ' isess{:} ' positive Corr'])
    set(gcf,'NextPlot','add');
    
    hfn(is) = figure; hold on
    suptitle(['session ' isess{:} ' negative Corr'])
    set(gcf,'NextPlot','add');
    
    for iu = unique(InfoTable(strcmp(InfoTable.Session,isess{:}),:).cluname,'stable')'
        
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
                        
                        
                        %% Categorize by temporal response or not
                        % Get unit data to filter by response quality
                        UnitData = UnitTable(ismember(UnitTable(:,1:12),T_pdc(1,[1:11 14]),'rows'),:);
                        if size(UnitData,1)>1, keyboard, end
                        
                        if isempty(UnitData) 
                            warning('  *** %s %s not found in DataTable! ',isess{:},iu{:})
                            continue
                        end
                        
                        if UnitData.Corr0 < 0
                            dpcat = 1; %negative
                        elseif isnan(UnitData.Corr0)
                            dpcat = 2; %insig
                        elseif UnitData.Corr0 >0
                            dpcat = 3; %positive
                        end
                        
                        %%
                        
                        % Calculate response for periodic stim
                        resp_pdc = AmalgamatedSpikes(ismember(InfoTable,T_pdc,'rows'),:);
                        plot_this = mean(binspikecounts(resp_pdc,binsize)/binsize*1000,1); %avg FR binned
                        
                        % Skip data points that likely come from junk clusters
                        if mean(plot_this)<5 || size(resp_pdc,1)<18
                            continue
                        end
                                                
                        
                        % Create empty matrices for population data if they
                        % don't exist already
                        if ~exist('Spikes_neg_j_0','var')
                            Spikes_neg_j_0=[];
                        end
                        if ~exist('Spikes_ins_j_0','var')
                            Spikes_ins_j_0=[];
                        end
                        if ~exist('Spikes_pos_j_0','var')
                            Spikes_pos_j_0=[];
                        end
                        
                        switch dpcat
                            case 1 %negative correlation
                                Spikes_neg_j_0 = [Spikes_neg_j_0; plot_this];
                                
                                % and plot unit
                                figure(hfn(is));
                                subplot(4,4,find(plot_order==str2double(extractAfter(strtok(iu{:},'_'),'ch'))));
                                hold on
                                set(gca,'XTick',[],'YTick',[])
                                plot(plot_this,'Color', [0 0 0],'LineWidth',4);
                                xlim([1 size(plot_this,2)])
                                title(iu{:})
                            
                            case 2 %insignificant correlation
                                Spikes_ins_j_0 = [Spikes_ins_j_0; plot_this];
                                
                            case 3 %positive correlations
                                Spikes_pos_j_0 = [Spikes_pos_j_0; plot_this];
                                
                                % and plot unit
                                figure(hfp(is));
                                subplot(4,4,find(plot_order==str2double(extractAfter(strtok(iu{:},'_'),'ch'))));
                                hold on
                                set(gca,'XTick',[],'YTick',[])
                                plot(plot_this,'Color', [0 0 0],'LineWidth',4);
                                xlim([1 size(plot_this,2)])
                                title(iu{:})
                        end
                        
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
                        plot_this = mean(binspikecounts(resp_jit,binsize)/binsize*1000,1); %avg FR binned
                        
                        % Skip data points that likely come from junk clusters
                        if mean(plot_this)<5 || size(resp_jit,1)<18
                            continue
                        end
                        
                        
                        % Create empty matrices for population data if they
                        % don't exist already
                        ij{:}(strfind(ij{:},'-'))=[];
                        if ~exist(sprintf('Spikes_neg_j_%s',ij{:}),'var')
                            eval(sprintf('Spikes_neg_j_%s = [];',ij{:}))
                        end
                        if ~exist(sprintf('Spikes_ins_j_%s',ij{:}),'var')
                            eval(sprintf('Spikes_ins_j_%s = [];',ij{:}))
                        end
                        if ~exist(sprintf('Spikes_pos_j_%s',ij{:}),'var')
                            eval(sprintf('Spikes_pos_j_%s = [];',ij{:}))
                        end
                        
                        switch dpcat
                            case 1 %negative correlation
                                eval(sprintf('Spikes_neg_j_%s = [Spikes_neg_j_%s; plot_this];',ij{:},ij{:}))
                                
                                % and plot unit
                                figure(hfn(is));
                                subplot(4,4,find(plot_order==str2double(extractAfter(strtok(iu{:},'_'),'ch'))));
                                hold on
                                set(gca,'XTick',[],'YTick',[])
                                plot(plot_this,...
                                    'Color',ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ),...
                                    'LineWidth',2);
                                
                            case 2 %insignificant correlation
                                eval(sprintf('Spikes_ins_j_%s = [Spikes_ins_j_%s; plot_this];',ij{:},ij{:}))
                                
                            case 3 %positive correlation
                                eval(sprintf('Spikes_pos_j_%s = [Spikes_pos_j_%s; plot_this];',ij{:},ij{:}))
                                
                                % and plot unit
                                figure(hfp(is));
                                subplot(4,4,find(plot_order==str2double(extractAfter(strtok(iu{:},'_'),'ch'))));
                                hold on
                                set(gca,'XTick',[],'YTick',[])
                                plot(plot_this,...
                                    'Color',ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ),...
                                    'LineWidth',2);
                        end
                        
                        
                    end
                end
            end
        end
    end %unit
end  %session


%% Plot grand average population PSTHs for each jitter condition

plot_mean = 1;

for iba = {'neg' 'ins' 'pos'}
    
    C = who('-regexp', ['Spikes_' iba{:}]);
    
    eval(sprintf('figure(hf_pop_%s);',iba{:}))
    hold on
    
    for ic = C'
        
        if isempty(eval(ic{:})), continue, end
        
        this_data=[]; outliers_rm=[];
        this_data = eval(ic{:});
        outliers_rm = this_data(abs(sum(zscore(this_data),2))<10,:);
        
        if ~plot_mean
            
            plot(repmat(1:(250/binsize),size(this_data,1),1)', this_data',...
                'Color', ALLcolors( strcmp( strtok(extractAfter(ic{:},[iba{:} '_j_']),'_'), strtok(plotOptions.colSelect,'_') ), : ),...
                'LineWidth',2)
        else
            
            if strcmp( strtok(extractAfter(ic{:},[iba{:} '_j_']),'_'), '0' )
                
%                 fill([1:(250/binsize) (250/binsize):-1:1],[mean(this_data,1) zeros(1,(250/binsize))],...
%                     [0 0 0],'FaceAlpha',0.4)
                
                fill([1:(250/binsize) (250/binsize):-1:1],[mean(this_data,1)+std(this_data,1) fliplr(mean(this_data,1)-std(this_data,1))],...
                    [0 0 0],'FaceAlpha',0.4)                
            else
                
                plot( mean(this_data,1) ,...
                    'Color', ALLcolors( strcmp( strtok(extractAfter(ic{:},[iba{:} '_j_']),'_'), strtok(plotOptions.colSelect,'_') ), : ),...
                    'LineWidth',6)
                
                [texty,textx] = max(mean(this_data,1));
                text(textx,texty,char(extractAfter(ic{:},[iba{:} '_j_'])))
            end
        end
    end
    
    
    set(gca,'XTick',[1:250/binsize]+0.5,...
        'XTickLabel',binsize:binsize:250)
    xlim([1 250/binsize])
    xlabel('ms during standard pd')
    ylabel('FR')
    
end



%% Save figures

savedir = fullfile(fn.standardPd,'PSTH_CategbyCorr0');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
 
% Population plots, below and above VS criterion
set(hf_pop_neg,'PaperOrientation','landscape');
print(hf_pop_neg,fullfile(savedir,'populationPSTH_NEG_Corr0'),'-dpdf','-bestfit');

set(hf_pop_ins,'PaperOrientation','landscape');
print(hf_pop_ins,fullfile(savedir,'populationPSTH_ins_Corr0'),'-dpdf','-bestfit');

set(hf_pop_pos,'PaperOrientation','landscape');
print(hf_pop_pos,fullfile(savedir,'populationPSTH_POS_Corr0_exTraces'),'-dpdf','-bestfit');


% Individual sessions/sites
for is = 1:numel(sessions)
    set(hfp(is),'PaperOrientation','landscape');
    print(hfp(is),fullfile(savedir,['sess' sessions{is} '_PSTHs_POS']),'-dpdf','-bestfit');
end

for is = 1:numel(sessions)
    set(hfn(is),'PaperOrientation','landscape');
    print(hfn(is),fullfile(savedir,['sess' sessions{is} '_PSTHs_NEG']),'-dpdf','-bestfit');
end








end






