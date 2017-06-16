function standardPd_populationPSTH(subject)

global fn 

%~~~~~~~~~~~~~~
binsize = 25;
%~~~~~~~~~~~~~~

fn = set_paths_directories;

% Load data
load(fullfile(fn.processed,'StandardPd_Spikes'));
InfoTable = readtable(fullfile(fn.processed,sprintf('StandardPd_StimInfo_%s',subject)));


%% Plot

% Set plot options 
[~, plotOptions] = setOptions;
plotOptions.colSelect = {'0_' '10_' '20_' '30_' '40_' '50_' '70_' '100_' '150_' '200_'};
ALLcolors = copper( numel(plotOptions.colSelect) );
set(0,'DefaultTextInterpreter','none')
set(0,'defaultaxesfontsize',plotOptions.labelSize)
set(0, 'DefaultFigureVisible', 'on')

% Get unique jitters
jitters = unique(InfoTable.jitter);
leg_str = cell(numel(jitters),1);
ils = 0;

% Plot
for iRate = unique(InfoTable.CenterRate)
    
    %% First plot Standard Pd response against other periods
    
    load(fullfile(fn.processed,'OtherPds_Spikes'));
    
    % Get indices of periodic stimuli/unit trials only
    pdc_idx = strcmp(InfoTable.jitter,'0');
    
    % Prepare plot
    PlotColors = hsv(size(OtherPdSpikes,3));
    leg_str2 = {'4' '2' '3' '5' '6' '7'};
    hfO = figure; hold on
    
    this_data = AmalgamatedSpikes(pdc_idx,:);
    this_data = binspikecounts(this_data,binsize)/binsize*1000;
    outliers_rm = this_data(abs(sum(zscore(this_data),2))<10,:);
    
    % Plot standard first
    plot_std = std(this_data,1);
    plot_this = mean(this_data,1); %avg FR binned
    
%     fill([1:size(plot_std,2) fliplr(1:size(plot_std,2))], [plot_this+plot_std plot_this-fliplr(plot_std)],...
%         [0 0 0],'FaceAlpha',0.5)
    ip(1) = plot(plot_this,'Color', [0 0 0],'LineWidth',7);
    
    for ipd = 1:size(OtherPdSpikes,3)
        
        this_data = OtherPdSpikes(pdc_idx,1:250,ipd);
        this_data = binspikecounts(this_data,binsize)/binsize*1000;
        outliers_rm = this_data(abs(sum(zscore(this_data),2))<10,:);
        
        plot_std = std(this_data,1);
        plot_this = mean(this_data,1);
        
%         fill([1:size(plot_std,2) fliplr(1:size(plot_std,2))], [plot_this+plot_std plot_this-fliplr(plot_std)],...
%             PlotColors(ipd,:),'FaceAlpha',0.3)
        ip(ipd+1)=plot(plot_this,'Color',PlotColors(ipd,:),'LineWidth',3);
        
    end
    legend(ip,leg_str2)
    xlim([1 size(plot_this,2)])
    set(gca,'XTick',[1:size(plot_this,2)]-0.5,...
        'XTickLabel',0:binsize:binsize*size(plot_this,2))
    xlabel('ms during standard pd')
    ylabel('average FR (sp/s)')
    title('Population average FR during periodic stim - middle "standard" period vs other periods')
    
    % Save figure
    set(hfO,'PaperOrientation','landscape');
    print(hfO,fullfile(fn.standardPd,'periodicPds_populationPSTH'),'-dpdf','-bestfit');
    
    
    %% Now plot jittered comparison
    hfJ = figure; hold on
    ipj=0;
    for ij = jitters'
        
        ipj=ipj+1;
        
        this_data = AmalgamatedSpikes(strcmp(InfoTable.jitter,ij),:);
        this_data = binspikecounts(this_data,binsize)/binsize*1000;
        outliers_rm = this_data(abs(sum(zscore(this_data),2))<10,:);
%         this_data = this_data(abs(sum(zscore(this_data),2))<10,:);
        
        plot_std = std(this_data,1);
        plot_this = mean(this_data,1); %avg FR binned
        
        % All together
        if strcmp(ij{:},'0')
            
            fill([1:size(plot_this,2) size(plot_this,2):-1:1],[plot_this zeros(1,size(plot_std,2))],...
                [0 0 0],'FaceAlpha',0.4)
        else
            %         fill([1:size(plot_std,2) fliplr(1:size(plot_std,2))], [plot_this+plot_std plot_this-fliplr(plot_std)],...
            %             ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ),...
            %             'FaceAlpha',0.3)
            ip(ipj)=plot(plot_this,'Color',...
                ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ),...
                'LineWidth',3);
            
            %         ip=plot(plot_this,...
            %             'Color', ALLcolors( strcmp(strtok(ij{:},'_'),strtok(plotOptions.colSelect,'_')), : ),...
            %             'LineWidth',3);
%             ip(ipj).LineWidth=10;
        end
        
        
        rV = jitter_LUT(iRate,ij{:});
        
        [texty,textx] = max(plot_this);
        text(textx,texty,[ij{:} ' *' num2str(rV(3)) '*  N=' num2str(sum(diff(InfoTable(strcmp(InfoTable.jitter,ij),:).trialN)<0)+1)])
        
        ils = ils+1;
        leg_str{ils,1} = [ij{:} ' *' num2str(rV(3)) '*  N=' num2str(sum(diff(InfoTable(strcmp(InfoTable.jitter,ij),:).trialN)<0)+1)];
    end
    
    xlim([1 size(plot_this,2)])
    set(gca,'XTick',[1:size(plot_this,2)]-0.5,...
        'XTickLabel',0:binsize:binsize*size(plot_this,2))
    xlabel('ms during standard pd')
    ylabel('average FR (sp/s)')
    title('Population average FR during standard 4 Hz period')
    legend(ip,leg_str,'Location','northeast','Interpreter','none')
    
    % Save figures
    set(hfJ,'PaperOrientation','landscape');
    print(hfJ,fullfile(fn.standardPd,'standardPd_populationPSTH'),'-dpdf','-bestfit');
    
    
    
    
end



end


