function ap_plot_FRresp_depths(subject,session,channel,clu,raster)

% TO DO:
%  mark drinking/active datapoints with different color
%  option for barplots or psychometric function
%  combine actual data across blocks
%  make even more modular, by measure type and plot type


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

if nargin<5
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end

% Find stimuli with more than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

blocks = unique([raster.block]);

% Make separate figures for each param except jitter and depth
for ib = blocks
    bk_raster = raster([raster.block]==ib);
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    % Calculate baseline rate for this block
    baselineFR = sum([bk_raster.x]<0) / (bk_raster(1).window_ms(1)/-1000) / sum(cellfun(@(x) ( numel(x) ), {bk_raster.tr_idx}));
    
for ip = 1:max(np)
    
    alldata = bk_raster(np==ip);
    [depths,~,ndpth] = unique([alldata.AMdepth]);
    
    
    hF = figure;
    hold on

    % Go through by increasing AMdepth
    for id = 1:max(ndpth)
        subplot(1,max(ndpth),id)
        hold on
        plot_subplot(alldata(ndpth==id),baselineFR)        
    end
    
    % Set axis limits to be same for all subplots
    hAllAxes = findobj(hF,'type','axes');
    if iscell(get(hAllAxes,'YLim'))
        yset(1) = min(cellfun(@min,get(hAllAxes,'YLim')));
        yset(2) = max(cellfun(@max,get(hAllAxes,'YLim')));
    else
        yset = get(hAllAxes(1),'YLim');
    end
    set(hAllAxes,'YLim',yset)
    
    % Add stimulus info in title (disregard jitter)
    suptitle(alldata(1).stim_str)
    
    % Get string for savename
    str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
    if numel(depths)==1
        str_dpth = num2str(depths*100);
    else
        str_dpth = 'ROVED';
    end
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    savename = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%i',...
        subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},ib);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    print(hF,'-depsc',fullfile(an_dir,savename))


end
end

% keyboard 

end

function plot_subplot(data,baselineFR)

% Set labels
if data(1).AMdepth==0
    jitter_labels = {repmat('unmodulated',1,numel(data))};
else
    jitter_labels = [data.jitter];
end

% Get behaving datapoints
drinking = find(strcmp({data.behaving},'D'));
behaving = find(strcmp({data.behaving},'A'));
if ~isempty(drinking)
    keyboard
end

jitters = str2double(strtok([data.jitter],'_'));
nj_max = sum(jitters==mode(jitters));

% Get FR and std of FR for error bars
FRvec = [data.nSpk];
FRstd = nan(1,numel(data)); FRtr=struct();
for ir = 1:numel(data)
    x = data(ir).x;
    y = data(ir).y;
    FR_t = nan(1,max(y));
    for it = 1:max(y)
        FR_t(it) = numel(x(x(y==it) > data(ir).AMonset & x(y==it) < data(ir).stimDur)) / ((data(ir).stimDur-data(ir).AMonset)/1000);
    end
    FRtr(ir).tr = FR_t;
    FRstd(ir) = std(FR_t,'omitnan');
end

% Set up vectors for creating bar plots
xbarvec = nan( 1,numel(unique(jitters))*nj_max);
ybarvec = nan( 1,numel(unique(jitters))*nj_max);
ebarvec = nan( 1,numel(unique(jitters))*nj_max);
xbarlab = cell(1,numel(unique(jitters))*nj_max);
barcols = copper(numel(unique(jitters)));

ii=1-nj_max;
ic=0;
for ij = unique(jitters)
    ii=ii+nj_max; ic=ic+1;
    nj = sum(jitters==ij);
    
    xbarvec(ii:ii+nj_max-1) = ij;
    
    if nj<nj_max
        ii1=ii+nj_max-nj-1;
        ybarvec(ii1:ii1+nj-1) = FRvec(jitters==ij);
        ebarvec(ii1:ii1+nj-1) = FRstd(jitters==ij);
        xbarlab(ii1:ii1+nj-1) = deal({jitter_labels{jitters==ij}});
        
        trFRs = {FRtr(jitters==ij).tr};
        ispc=0;
        for isp=ii1:ii1+nj-1
            ispc=ispc+1;
            plot(repmat(isp,size(trFRs{ispc})), trFRs{ispc}, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        end
        
    elseif nj==nj_max
        ybarvec(ii:ii+nj_max-1) = FRvec(jitters==ij);
        ebarvec(ii:ii+nj_max-1) = FRstd(jitters==ij);
        xbarlab(ii:ii+nj_max-1) = deal({jitter_labels{jitters==ij}});
        
        trFRs = {FRtr(jitters==ij).tr};
        ispc=0;
        for isp=ii:ii+nj_max-1
            ispc=ispc+1;
            plot(repmat(isp,size(trFRs{ispc})), trFRs{ispc}, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        end
        
    end
    
    % Plot this jitter group
    b=bar(ii:ii+nj_max-1, ybarvec(ii:ii+nj_max-1), 'BaseValue',baselineFR, 'FaceColor', barcols(ic,:));
    if ij==0
        b.EdgeColor = 'blue';
        b.LineWidth = 2;
    end
    errorbar(ii:ii+nj_max-1, ybarvec(ii:ii+nj_max-1), ebarvec(ii:ii+nj_max-1),...
        'LineStyle','none', 'Color',[0 0 0] , 'LineWidth',1 )
    
end

mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(data(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    ylabel('avg FR during AM (Hz)')
end



end

