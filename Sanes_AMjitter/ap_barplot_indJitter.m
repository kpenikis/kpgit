function ap_barplot_indJitter(subject,session,channel,clu,METRIC,raster)
% Plots various response measures for each stimulus, as defined in input
% variable PLOT_TYPE. Current options are FR and FF, as bar plots only.
% CV = std / mean;
% FF = var / mean;

% TO DO:
%  mark drinking/active datapoints with different color
%  option for barplots or psychometric function
%  combine actual data across blocks
%  make even more modular, by measure type and plot type


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath('analysis_modules','helpers')

if nargin<6 || ~exist('raster','var')
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end



% Extra param to set for Fano Factor analysis
binsize = 250;

% Find stimuli with more than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


% Find blocks and designate which ones to combine
blocks = unique([raster.block]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
group_blocks = [89 90];  %only 2 at a time for now
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ic = size(group_blocks,1)
    blocks(blocks==group_blocks(ic,2)) = [];
end

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
    
    % Calculate baseline rate for this block
    baselineFR = sum([bk_raster.x]<0) / (bk_raster(1).window_ms(1)/-1000) / sum(cellfun(@(x) ( numel(x) ), {bk_raster.tr_idx}));
    
    
    for ip = 1:max(np)
        
        param_raster = bk_raster(np==ip);
        
        [depths,~,ndpth] = unique([param_raster.AMdepth]);
        
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
                case 'VS'
                    data_mean = calc_VS(stim,subject);
                    data_std = nan(length(data_mean),2); data_trs = nan;
                case 'RS'
                    data_mean = calc_RS(stim,subject);
                    data_std = nan(length(data_mean),2); data_trs = nan;
            end
            pause(0.2)
            subplot_bar(stim,data_mean,data_std,data_trs,baselineFR,METRIC);
            pause(0.2)
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
            str_dpth = 'ROVED';
        end
        
        % Add stimulus info in title
        title_str = sprintf('ch %i clu%i\nnoise: %s - %s Hz  |  %sdB\nAM %s Hz  |  %s dpth  |  blk%s',...
            channel,clu,str_pars{1},str_pars{2},str_pars{3},str_pars{4},str_dpth,bk_str);
        suptitle(title_str);
        
        % Save figure
        an_dir = fullfile(savedir,subject,'^an_plots',session);
        if ~exist(an_dir,'dir')
            mkdir(an_dir)
        end
        
        switch METRIC
            case 'FR'
                savename = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'FF'
                savename = sprintf('%s_%s_FanoFactor_bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'VS'
                savename = sprintf('%s_%s_VS_bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
            case 'RS'
                savename = sprintf('%s_%s_Rayleigh_bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                    subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        end
        
        print(hF,'-depsc',fullfile(an_dir,savename))
        
        
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

% Draw poisson lie for FF
if any(strcmp(METRIC,{'VS' 'FF'}))
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
        b=bar(ii, ybarvec(ii), 'BaseValue', baselineFR, 'FaceColor', barcols(ic,:));
        
        % Plot FR from individual trials
        trFRs = data_trs(ii).tr;
        plot(repmat(ii,size(trFRs)), trFRs, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        
    else
        b=bar(ii, ybarvec(ii), 'FaceColor', barcols(ic,:));
        
    end
    
    if drinking(ii)==1
        b.FaceAlpha = 0.4;
        b.EdgeColor = barcols(ic,:);
        b.LineWidth = 2;
    elseif behaving(ii)==1
        b.FaceAlpha = 0.7;
        b.EdgeColor = 'r';
        b.LineWidth = 1;
    end
    if jitters(ii)==0
        b.EdgeColor = 'blue';
        b.LineWidth = 2;
    end
    
    
    errorbar(ii, ybarvec(ii), ebarvec(1,ii), ebarvec(2,ii),...
        'LineStyle','none', 'Color',[0 0 0] , 'LineWidth',1 )

end

mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');

% bk_str = num2str(unique([stim.block]));
% blocks = unique([stim.block]);
% for ib = 1:numel(unique([param_raster.block]))
%     bk_str = [bk_str ];
% end

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(stim(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    switch METRIC
        case 'FR'
            ylabel('avg FR during AM')
        case 'FF'
            ylabel('avg FF during AM')
        case 'VS'
            ylabel('avg Vector Strength')
        case 'RS'
            ylabel('avg Rayleigh Statistic')
    end
else
    set(gca,'YTickLabel',[])
end


end % subplot_bar




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


