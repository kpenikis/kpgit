function ap_analysis_plots_v2(subject,session,channel,clu,METRIC,PLOT_TYPE,raster)
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

if nargin<7 || ~exist('raster','var')
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end


% for j = 1:(length(varargin)/2)
%     varargin{(2*j)-1 }
%     varargin{2*j} 
% end


% Find stimuli with more than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

% Extra params to set, for Fano Factor analysis
binsize = 250;


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
        if id==1, continue, end
        % Collapse data across blocks, if params repeated
        stim = collapse_identical_blocks(param_raster(ndpth==id));
        
        % Plot subplot
        subplot(1,max(ndpth),id)
        hold on
        switch METRIC
            case 'FR'
                [data_mean,data_std,data_trs] = calc_FR(stim);
            case 'FF'
                [data_mean,data_std] = calc_FF(stim,binsize);
                data_trs = nan;
            case 'FFpd'
                [data_mean] = calc_FF_periods(stim,subject);
                data_std = nan; data_trs = nan;
        end
        
        switch PLOT_TYPE
            case 'bar'
                subplot_bar(stim,data_mean,data_std,data_trs,baselineFR);
            case 'matrix'
                imagesc(data_mean)
                colorbar
        end
        
    end
    
    % Set axis limits to be same for all subplots
    hAllAxes = findobj(hF,'type','axes');
    if iscell(get(hAllAxes,'YLim'))
        ymax = max(cellfun(@max,get(hAllAxes,'YLim')));
    else
        ymax = max(get(hAllAxes(1),'YLim'));
    end
    set(hAllAxes,'YLim',[0 ymax])
    
    % Add stimulus info in title (disregard jitter)
    suptitle(param_raster(1).stim_str)                       %%%%% change this too
    
    % Get string for savename
    str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
    if numel(depths)==1
        str_dpth = num2str(depths*100);
    else
        str_dpth = 'ROVED';
    end
    
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    
    switch PLOT_TYPE
        case 'FR'
            savename = sprintf('%s_%s_FRresp_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
        case 'FF'
            savename = sprintf('%s_%s_FanoFactor_bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
                subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);
    end
    
%     print(hF,'-depsc',fullfile(an_dir,savename))


end
end

% keyboard 

end




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function stim = collapse_identical_blocks(raster)

% Combine actual data across blocks
drinking = strcmp({raster.behaving},'D');
[unq,~,iu] = unique([[raster.fileIDs]; drinking]','rows');
for iiu = 1:max(iu)
    
    rs = find([raster.fileIDs]==unq(iiu,1) & drinking==unq(iiu,2));
    
    x=[]; y=0; bk=[]; 
    for irs = rs
        x = [x raster(irs).x];
        y = [y raster(irs).y + max(y)*ones(size(raster(irs).y))];
        bk = [bk raster(irs).block];
    end
    y(1)=[];
    
    stim(iiu) = raster(rs(1));
    stim(iiu).block = bk;
    stim(iiu).x = x; 
    stim(iiu).y = y;
    
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function [FRvec,FRstd,FRtr] = calc_FR(stim)

% Get FR and std of FR for error bars
FRvec = nan(1,numel(stim)); FRstd = nan(1,numel(stim)); FRtr=struct();
for ir = 1:numel(stim)
    x = stim(ir).x;
    y = stim(ir).y;
    FR_t = nan(1,max(y));
    for it = 1:max(y)
        FR_t(it) = numel(x(x(y==it) > stim(ir).AMonset & x(y==it) < stim(ir).stimDur)) / ((stim(ir).stimDur-stim(ir).AMonset)/1000);
    end
    FRtr(ir).tr = FR_t;
    FRvec(ir) = mean(FR_t,'omitnan');
    FRstd(ir) = std(FR_t,'omitnan');
end

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [FFmean,FFstd] = calc_FF(stim,binsize)

% Find minimum number of trials for each stimulus
ntc=nan(numel(stim),1);
for is = 1:numel(stim)
    ntc(is) = max(stim(is).y);
end
min_nt = min(ntc);
% min_nt = 5;

iterations = 200;  rng('shuffle');

% Set up empty vectors and get FF data
FF     = nan(numel(stim),iterations);

for ii = 1:iterations
for is = 1:numel(stim)
    
    data = stim(is);
    
    if binsize==0
        tVec = [data.AMonset data.stimDur];
    else
        tVec = data.AMonset : binsize : data.stimDur;
    end
    
    trs = 1:max(data.y); 
    trs(randi(max(data.y),[max(trs)-min_nt 1])) = []; %remove random trials
    sp_hist = nan(min_nt,length(tVec)-1);
    
    ir = 0;
    for it = trs
        ir=ir+1;
        
        sp = data.x(data.y==it);
        sp_hist(ir,:) = histcounts(sp,tVec);
        
    end
    
    % Calculate means and variances for time bins
    tr_var  = var(sp_hist,1);
    tr_mean = mean(sp_hist,1);
    
    % Calculate Fano Factor for this stimulus
    FF(is,ii) = mean(tr_var(tr_mean~=0) ./ tr_mean(tr_mean~=0));

end
end

FFmean = mean(FF,2);
FFstd = std(FF,1,2)/sqrt(iterations);

end


function [FFmean] = calc_FF_periods(stim,subject)

% Set stimulus file directory
blocks = stim.block;
stimdir = fullfile('/Users/kpenikis/Documents/SanesLab/Data/raw_data',subject,sprintf('Block-%i_Stim',blocks(1)));
rV = load(fullfile(stimdir,stim(1).stimfn));

% Find minimum number of trials for each stimulus
ntc=nan(numel(stim),1);
for is = 1:numel(stim)
    ntc(is) = max(stim(is).y);
end
min_nt = min(ntc);
% min_nt = 5;

iterations = 200;  rng('shuffle');

% Set up empty vectors and get FF data
FF = nan( numel(stim), length(rV.buffer)-1, iterations );

for is = 1:numel(stim)
    
    data = stim(is);
    
    % Get vectors of rates for this stimulus
    rateVec = load(fullfile(stimdir,data.stimfn));
    rateVec = rateVec.buffer;
    tVec = round(data.AMonset + cumsum([0 0.75*(1000/rateVec(2)) 1000./rateVec(3:end)]));
    
for ii = 1:iterations
    
    trs = 1:max(data.y); 
    trs(randi(max(data.y),[max(trs)-min_nt 1])) = []; %remove random trials
    sp_hist = nan(min_nt,length(tVec)-1);
    
    ir = 0;
    for it = trs
        ir=ir+1;
        
        sp = data.x(data.y==it);
        sp_hist(ir,:) = histcounts(sp,tVec);
        
    end
    
    % Calculate means and variances for time bins
    tr_var  = var(sp_hist,1);
    tr_mean = mean(sp_hist,1);
    
    % Calculate Fano Factor for this stimulus
    FF(is,:,ii) = tr_var(tr_mean~=0) ./ tr_mean(tr_mean~=0);

end
end

FFmean = mean(FF,3);
FFstd  = std(FF,1,3)/sqrt(iterations);

end




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%        SUBPLOTS 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function subplot_bar(stim,data_mean,data_std,data_trs,baselineFR)


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

% Set up vectors for creating bar plots
xbarvec = nan( 1,numel(stim));
ybarvec = nan( 1,numel(stim));
ebarvec = nan( 1,numel(stim));
xbarlab = cell(1,numel(stim));
barcols = copper(numel(stim)); %numel(unique(jitters)) because of drinking

ic=0;
for ii = 1:numel(stim)
    ic=ic+1;
    
    % Get data to plot
    xbarvec(ii) = jitters(ii);
    ybarvec(ii) = data_mean(ii);
    ebarvec(ii) = data_std(ii);
    xbarlab(ii) = deal({jitter_labels{ii}});
    
    % Plot some things only if this is a FR plot
    if ~isnan(data_trs) 
        b=bar(ii, ybarvec(ii), 'BaseValue', baselineFR, 'FaceColor', barcols(ic,:));
        
        % Plot FR from individual trials
        trFRs = data_trs(ii).tr;
        plot(repmat(ii,size(trFRs)), trFRs, 'o', 'Color', barcols(ic,:),'MarkerSize',8)
        
    else
        b=bar(ii, ybarvec(ii), 'FaceColor', barcols(ic,:));
    end
    
    if jitters(ii)==0
        b.EdgeColor = 'blue';
        b.LineWidth = 2;
    end
    if drinking(ii)==1
        b.FaceAlpha = 0.2;
        b.EdgeColor = barcols(ic,:);
        b.LineWidth = 2;
    elseif behaving(ii)==1
        b.FaceAlpha = 0.5;
        b.EdgeColor = 'r';
        b.LineWidth = 1;
    end
    errorbar(ii, ybarvec(ii), ebarvec(ii),...
        'LineStyle','none', 'Color',[0 0 0] , 'LineWidth',1 )
    
end

mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(stim(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    ylabel('avg FR during AM (Hz)')
end


end % function




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plot_subplot_FF(raster,binsize)

% Set labels
if raster(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(raster));
else
    jitter_labels = [raster.jitter];
end
jitters = str2double(strtok([raster.jitter],'_'));
nj_max = sum(jitters==mode(jitters));

% Get behaving datapoints
drinking = find(strcmp({raster.behaving},'D'));
behaving = find(strcmp({raster.behaving},'A'));
if ~isempty(drinking)
%     keyboard
end

% Merge raster data across blocks
%   have to change substantially -- get FF data by stepping through
%   jitters, as opposed to entries in data structure
% 1) collapse data from blocks that are identical in all other params,
% including Beh
% 2) color according to Beh

% Set up empty vectors and get FF data
FF     = nan(numel(raster),1);
FFtime = nan(numel(raster),500);

% Find minimum number of trials for each stimulus
ntc=nan(numel(raster),1);
for ii = 1:numel(raster)
    ntc(ii) = max(raster(ii).y);
end
ent = min(ntc);

for ii = 1:numel(raster)
    
    data = raster(ii);
    
    if binsize==0
        tVec = [data.AMonset data.stimDur];
    else
        tVec = data.AMonset : binsize : data.stimDur;
    end
    
    trs = 1:max(data.y); 
    trs(randi(max(data.y),[max(trs)-ent 1])) = []; %remove random trials
    sp_hist = nan(ent,length(tVec)-1);
    
    ir = 0;
    for it = trs
        ir=ir+1;
        
        sp = data.x(data.y==it);
        sp_hist(ir,:) = histcounts(sp,tVec);
        
    end
    
    % Calculate means and variances for time bins
    tr_var  = var(sp_hist,1);
    tr_mean = mean(sp_hist,1);
    
    % Calculate Fano Factor for this stimulus
    FF(ii) = mean(tr_var(tr_mean~=0) ./ tr_mean(tr_mean~=0));
    FFtime(ii,1:length(tr_mean)) = tr_var ./ tr_mean;

end


%%%%%%%%%%%%%%%%%%%%%%%
% Now prepare the plot 

% Set up vectors for creating bar plots
xbarvec = nan( 1,numel(unique(jitters))*nj_max);
ybarvec = nan( 1,numel(unique(jitters))*nj_max);
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
        ybarvec(ii1:ii1+nj-1) = FF(jitters==ij);
        xbarlab(ii1:ii1+nj-1) = deal({jitter_labels{jitters==ij}});
        
    elseif nj==nj_max
        ybarvec(ii:ii+nj_max-1) = FF(jitters==ij);
        xbarlab(ii:ii+nj_max-1) = deal({jitter_labels{jitters==ij}});
        
    end
    
    % Plot this jitter group
    b=bar(ii:ii+nj_max-1, ybarvec(ii:ii+nj_max-1), 'FaceColor', barcols(ic,:));
    if ij==0
        b.EdgeColor = 'blue';
        b.LineWidth = 2;
    end
    hold on

end

mean_j0FF = mean(ybarvec(xbarvec==0),'omitnan');

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
plot([0 length(xbarvec)+1],[mean_j0FF mean_j0FF],'--b')
plot([0 length(xbarvec)+1],[1 1],'--k')
set(gca,'xlim',[0 length(xbarvec)+1])
title([num2str(data(1).AMdepth*100) '% depth' ])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    ylabel('avg Fano Factor across time during AM')
end


end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







