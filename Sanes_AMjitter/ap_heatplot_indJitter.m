function ap_heatplot_indJitter(subject,session,channel,clu,METRIC,raster)
% Plots various response measures for each stimulus, as defined in input
% variable PLOT_TYPE. Current options are FR and FF, as bar plots only.
% CV = std / mean;
% FF = var / mean;

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath('analysis_modules','helpers')

if nargin<6 || ~exist('raster','var')
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end


% Remove stimuli with fewer than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

% Extra params to set, for Fano Factor analysis
binsize = 250;


% Find blocks and designate which ones to combine
blocks = unique([raster.block]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
group_blocks = [89 90];  %only 2 at a time for now      %make a function containing list of blocks to merge
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
    
    plot_data=[];
    % Go through by increasing AMdepth
    for id = 1:max(ndpth)
        
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
    switch METRIC
        case 'RS'
            cmax = 25;
        case 'VS'
            cmax = 1;
        case 'FF'
            cmax = 2;
    end
    hAllAxes = findobj(hF,'type','axes');
    set(hAllAxes,'CLim',[0 cmax])
    
    
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
%         case 'FR'
%             keyboard
%             savename = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
%                 subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'FF'
            savename = sprintf('%s_%s_FFpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'VS'
            savename = sprintf('%s_%s_VSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'RS'
            savename = sprintf('%s_%s_RSpd_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
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

% Plot data
imagesc(data);
colormap('gray')

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
xlim([0.5 size(data,2)+0.5])
ylim([0.5 size(data,1)+0.5])
xlabel('sAM period (chronological)')
title([num2str(stim(1).AMdepth*100) '% depth'])
% Color bar
h=colorbar; ylabel(h,METRIC,'FontSize',14)



end % subplot_mtx




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


