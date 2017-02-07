function ap_heatplot_indJitter(subject,session,channel,clu,METRIC)
% Plots various response measures for each stimulus, as defined in input
% variable METRIC. Calculates response for each AM period separately, and
% plots as a heat plot.


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers')


% Load raster struct
raster = get_raster(subject,session,channel,clu);

% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);


% Get blocks for loop and designate which ones to combine
[blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));


% Extra params to set, for Fano Factor analysis
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
    
    % Calculate baseline rate for this block
    baselineFR = sum([bk_raster.x]<0) / (bk_raster(1).window_ms(1)/-1000) / sum(cellfun(@(x) ( numel(x) ), {bk_raster.tr_idx}));
    
for ip = 1:max(np)
    
    param_raster = bk_raster(np==ip);
    [depths,~,ndpth] = unique([param_raster.AMdepth]);
    
    
    hF = figure;
    hold on
    
    plot_data=[];
    % Go through by increasing AMdepth
    for id = max(ndpth):-1:1
        
        % Collapse data across blocks, if params repeated
        stim = collapse_blocks(param_raster(ndpth==id));
        
        % Plot subplot
        nsp = [ (1+sum(ndpth<id)) : (sum(ndpth<id)+sum(ndpth==id)) ];
        subplot(length(ndpth), 1, nsp ,'align')
        hold on
        
        switch METRIC
            case 'FF'
                [data] = calc_FF_periods(stim,subject);
                data_std = nan; data_trs = nan;
            case 'VS'
                [~,data] = calc_VS(stim,subject);
%                 data_std = nan(length(data_mean),2); data_trs = nan;
            case 'RS'
                [~,data] = calc_RS(stim,subject);
%                 data_std = nan(length(data_mean),2); data_trs = nan;
        end
        
        pause(0.2)
        subplot_mtx(stim,data,METRIC);
        pause(0.2)
        
    end
    
    
    % Set COLOR axis limits to be same for all subplots
%     hAllAxes = findobj(hF,'type','axes');
%     switch METRIC
%         case 'RS'
%             cmax = 25;
%         case 'VS'
%             cmax = 1;
%             
%         case 'FF'        %automated axis limit
%             % Get current max val
%             if iscell(get(hAllAxes,'CLim'))
%                 cmax1 = max(cellfun(@max,get(hAllAxes,'CLim')));
%             else
%                 cmax1 = max(get(hAllAxes(1),'CLim'));
%             end
%             % Snap to closest integer 
%             if cmax1<4.5
%                 cmax = ceil(cmax1+0.3);
%             elseif cmax1>=4.5
%                 cmax = 6;
%                 keyboard
%             end
%     end
%     set(hAllAxes,'CLim',[0 cmax])
    
    
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
    
    
    
    %% Save figure
    
    switch METRIC
        case 'FF'
            savefolder = 'FF-periods';
            savename   = sprintf('%s_%s_FFpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'VS'
            savefolder = 'VS-periods';
            savename   = sprintf('%s_%s_VSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'RS'
            savefolder = 'RS-periods';
            savename   = sprintf('%s_%s_RSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
    end
    
    datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
    an_dir = fullfile(savedir,subject,'^an_plots',session,savefolder);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    
    set(gcf,'PaperOrientation','landscape');
    print(hF,'-dpdf',fullfile(an_dir,savename),'-bestfit')
%     print(hF,'-depsc',fullfile(an_dir,savename))


end
end

% keyboard 

end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        SUBPLOTS 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function subplot_mtx(stim,data,METRIC)


% Set color categories depending on metric
C=nan(size(data));
switch METRIC
    case 'FF'
        b1 = 0.4;
        b2 = 0.8;
        b3 = 1.2;
        b4 = 1.6;
    case 'VS'
        b1 = 0.15;
        b2 = 0.3;
        b3 = 0.45;
        b4 = 0.6;
    case 'RS'
        b1 = 6.9;
        b2 = 13.8;
        b3 = 20.7;
        b4 = 34.5;
end
C(data<b1)             = 1;
C(data>=b1 & data<b2)  = 2;
C(data>=b2 & data<b3)  = 3;
C(data>=b3 & data<b4)  = 4;
C(data>=b4)            = 5;

% Create color map
% map = [ 0.4  0.8  1     ;...
%         0.48 0.7  0.8   ;...
%         0.58 0.58 0.58  ;...
%         0.63 0.5  0.4   ;...
%         0.7  0.4 0.2   ];
map = [ 0.7  0.4  0.2  ;...
        0.85 0.7  0.6  ;...
        0.9  0.9  0.9  ;...
        0.55 0.55 0.55 ;...
        0.2  0.2  0.2  ];
colormap(flipud(map))

% Plot data
image(C)


% Create labels
if stim(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(stim));
else
    jitter_labels = [stim.jitter];
end
ylab = cell(1,numel(stim));
for ii = 1:numel(stim)
    
    ylab(ii) = deal({jitter_labels{ii}});
    
    if strcmp({stim(ii).behaving},'D')
        ylab{ii} = sprintf('%s: D',ylab{ii});
    end
    if strcmp({stim(ii).behaving},'A')
        ylab{ii} = sprintf('%s: A',ylab{ii});
    end
    
end
% Set labels
set(gca,'YTick',1:numel(stim),'YTickLabel',ylab,'TickLabelInterpreter','none')

% Other figure properties
title([num2str(stim(1).AMdepth*100) '% depth'])
xlim([0.5 size(data,2)+0.5])
ylim([0.5 size(data,1)+0.5])
if numel(findobj(gcf,'type','axes'))==1
    xlabel('sAM period (chronological)')
else
    set(gca,'XTick',[])
end
% Color bar
h=colorbar; 
h.Ticks = [1 2 3 4 5];
h.TickLabels = {'0' num2str(b1) num2str(b2) num2str(b3) num2str(b4)};
ylabel(h,METRIC,'FontSize',14)

% Plot horizontal lines between conditions
ylms = get(gca,'YLim');
xlms = get(gca,'XLim');
plot(repmat(xlms,diff(ylms)-1,1)',[(ylms(1)+1):(ylms(2)-1); (ylms(1)+1):(ylms(2)-1)],'k')



end % subplot_mtx




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


