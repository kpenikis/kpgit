function ap_barplot_indJitter(subject,session,channel,clu,METRIC)
% Calculates various response measures for each stimulus, as defined in input
% variable METRIC, and plots the results as bar plots.

%  (FF = var / mean)


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers')

% Load raster struct
raster = get_raster(subject,session,channel,clu);

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


% Get blocks for loop and designate which ones to combine
[blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));

% Extra param to set for Fano Factor analysis
binsize = 250;

% Make separate figures for each param except jitter and depth
for ib = blocks
    
    % Get data for these blocks
    bk_raster = raster([raster.block]==ib);
    
    % Add data from blocks set to combine, if needed
    if ~isempty(group_blocks(:,1)==ib)
        bk_raster = [bk_raster raster([raster.block]== group_blocks(group_blocks(:,1)==ib,2) )];
        bk_str = [num2str(ib) num2str(group_blocks(group_blocks(:,1)==ib,2))];
    else
        bk_str = num2str(ib);
    end
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    % Calculate baseline rate for this condition
% %     baselineFR = calc_baselineFR(bk_raster);
    
    
    for ip = 1:max(np)
        
        param_raster = bk_raster(np==ip);
        
        [depths,~,ndpth] = unique([param_raster.AMdepth]);
        
        % Calculate baseline rate for this condition
        baselineFR = calc_baselineFR(param_raster);
        
        hF = figure;
        hold on
        
        % Go through by increasing AMdepth
        for id = 1:max(ndpth)
            
            % Collapse data across blocks, if params repeated
            stim = collapse_blocks(param_raster(ndpth==id));
            
            % Plot subplot
            nsp = [ (1+sum(ndpth<id)) : (sum(ndpth<id)+sum(ndpth==id)) ];
            subplot(1,length(ndpth), nsp ,'align')
            hold on
            switch METRIC
                case 'FR'
                    [data_mean,data_std,data_trs] = calc_FR(stim);
                case 'FF'
                    [data_mean,data_std] = calc_FF(stim,binsize);
                    data_trs = nan;
                case 'FFpd-avg'
                    [data_trs] = calc_FF_periods(stim,subject);
                    data_mean = mean(data_trs,2);
                    data_std  = [std(data_trs,0,2) std(data_trs,0,2)]./2;
                case 'VS'
                    data_mean = calc_VS(stim,subject);
                    data_std = nan(length(data_mean),2); data_trs = nan;
                case 'RS'
                    data_mean = calc_RS(stim,subject);
                    data_std = nan(length(data_mean),2); data_trs = nan;
                case 'standardFR'
                    [data_mean,data_std,data_trs] = calc_standardFR(stim,subject);
            end
            
            % Plot data in subplot
            pause(0.1)
            subplot_bar(stim,data_mean,data_std,data_trs,baselineFR,METRIC);
            pause(0.4)
            
        end
        
        % Set axis limits to be same for all subplots
        hAllAxes = findobj(hF,'type','axes');
        if iscell(get(hAllAxes,'YLim'))
            ymax = max(cellfun(@max,get(hAllAxes,'YLim')));
        else
            ymax = max(get(hAllAxes(1),'YLim'));
        end
        set(hAllAxes,'YLim',[0 ymax])
        
        % Get string for savename
        str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
        if numel(depths)==1
            str_dpth = num2str(depths*100);
        else
            str_dpth = 'VARIED';
        end
        
        % Add stimulus info in title
        title_str = sprintf('ch %i clu%i\nnoise: %s - %s Hz  |  %sdB\nAM %s Hz  |  %s dpth  |  blk%s',...
            channel,clu,str_pars{1},str_pars{2},str_pars{3},str_pars{4},str_dpth,bk_str);
        suptitle(title_str);
        
        
        % Save figure
        
        switch METRIC
            case 'FR'
                savefolder = 'FRavg';
                savename   = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FF'
                savefolder = ['FFavg-' num2str(binsize)];
                savename   = sprintf('%s_%s_FF-bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FFpd-avg'
                savefolder = 'FFavg-periods';
                savename   = sprintf('%s_%s_FF-pdAVG_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'VS'
                savefolder = 'VSavg';
                savename   = sprintf('%s_%s_VS_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'RS'
                savefolder = 'RSavg';
                savename   = sprintf('%s_%s_Rayleigh_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'standardFR'
                savefolder = 'standardFR';
                savename   = sprintf('%s_%s_standardFR_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        end
        
        datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
        an_dir = fullfile(datadir,subject,'^an_plots',session,savefolder);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        set(gcf,'PaperOrientation','landscape');
        print(hF,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%         print(hF,'-depsc',fullfile(an_dir,savename))
        
pause(1)
    end
end

% keyboard 

end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        SUBPLOTS 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function subplot_bar(stim,data_mean,data_std,data_trs,baselineFR,METRIC)


% Set labels
if stim(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(stim));
else
    jitter_labels = [stim.jitter];
end
jitters = str2double(strtok([stim.jitter],'_'));

% Get behaving datapoints
drinking = strcmp({stim.behaving},'D');
behaving = strcmp({stim.behaving},'A');

% Recalculate baseline for different conditions
if sum([drinking,behaving])>0
    baseFR = nan(size(stim));
    
    st_d = stim(drinking   &  ~behaving );
    st_p = stim( ~drinking &  ~behaving );
    st_a = stim( ~drinking &  behaving  );
    
    baseFR(drinking)              = calc_baselineFR(st_d);
    baseFR(behaving)              = calc_baselineFR(st_a);
    baseFR(~behaving & ~drinking) = calc_baselineFR(st_p);
    
else
    baseFR = baselineFR*ones(size(stim));
end

% Draw poisson lie for FF
if any(strcmp(METRIC,{'VS' 'FF' 'FFpd-avg'}))
    plot([0 numel(stim)+1],[1 1],'-','Color',[0.5 0.5 0.5])
elseif strcmp(METRIC,'RS')
    plot([0 numel(stim)+1],[13.8 13.8],'-','Color',[0.5 0.5 0.5])
end

% Set up vectors for creating bar plots
xbarvec = nan( 1,numel(stim));
ybarvec = nan( 1,numel(stim));
ebarvec = nan( 2,numel(stim));
xbarlab = cell(1,numel(stim));
barcols = copper(numel(stim)); %numel(unique(jitters)) because of drinking

ic=0;
for ii = 1:numel(stim)
    ic=ic+1;
    
    % Get data to plot
    xbarvec(ii) = jitters(ii);
    ybarvec(ii) = data_mean(ii);
    ebarvec(:,ii) = data_std(ii,:)';
    xbarlab(ii) = deal({jitter_labels{ii}});
    
    % Plot some things only if this is a FR plot
    if strcmp(METRIC,'FR')
        b=fill([-0.4 -0.4 0.4 0.4]+ii, [baseFR(ii) ybarvec(ii) ybarvec(ii) baseFR(ii)], barcols(ic,:));
    else
        b=bar(ii, ybarvec(ii), 'FaceColor', barcols(ic,:),'EdgeColor',barcols(ic,:));
    end
    
    switch METRIC
        % Plot FR from individual trials
        case 'FR'
            trFRs = data_trs(ii).tr;
            plot(repmat(ii,size(trFRs)), trFRs, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        case 'FFpd-avg'
%             plot(repmat(ii,size(data_trs,2)), data_trs(ii,:), 'o', 'Color', barcols(ic,:),'MarkerSize',8)
    end
    
    % Special formatting by condition
    b.EdgeColor = barcols(ic,:);
    b.LineWidth = 1;
    if drinking(ii)==1
        b.FaceAlpha = 0.7;
    elseif behaving(ii)==1
        b.FaceAlpha = 0.4;
        b.EdgeColor = 'r';
    end
    if jitters(ii)==0
        b.EdgeColor = 'blue';
    end
    if stim(1).AMdepth==0
        b.EdgeColor = [0 0.4 0];
        b.LineWidth = 2;
    end
    
    % Add error bars
    errorbar(ii, ybarvec(ii), ebarvec(1,ii), ebarvec(2,ii),...
        'LineStyle','none', 'Color',[0 0 0], 'LineWidth',1 )

end

% Plot mean baseline FR and jitter=0 FR
if strcmp(METRIC,'FR')
    mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');
    plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
    plot([0 length(xbarvec)+1], [mean(baseFR) mean(baseFR)],'Color',[0.4 0.4 0.4],'LineWidth',0.2)
end

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(stim(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    switch METRIC
        case 'FR'
            ylabel('avg FR during AM')
        case 'FF'
            ylabel('avg FF during AM')
        case 'FFpd-avg'
            ylabel('avg FF during AM, by periods')
        case 'VS'
            ylabel('avg Vector Strength')
        case 'RS'
            ylabel('avg Rayleigh Statistic')
        case 'standardFR'
            ylabel('avg FR during standard period')
    end
else
    set(gca,'YTickLabel',[])
end


end % subplot_bar




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


