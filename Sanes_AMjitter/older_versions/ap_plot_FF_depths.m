function ap_plot_FF_depths(subject,session,channel,clu,raster)
% Plots a folded histogram of spiking during the AM portion of a stimulus.
% Note that periods are not all of equal length so n trials tapers off.
% CV = std / mean;
% FF = var / mean;

% TO DO:
%  mark drinking/active datapoints with different color
%  replace barplots with psychometric functions for roved depth, maybe all


set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

if nargin<5
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end

% Find stimuli with more than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

% Extra params to set, for Fano Factor analysis
binsize = 250;

% Find blocks and designate which ones to combine
blocks = unique([raster.block]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
combine_blocks = [89 90];  %only 2 at a time for now
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for ic = size(combine_blocks,1)
    blocks(blocks==combine_blocks(ic,2)) = [];
end



% Make separate figures for each param except jitter and depth
for ib = blocks
    
    % Get data for these blocks
    bk_raster = raster([raster.block]==ib);
    
    % Add data from blocks set to combine, if needed
    if ~isempty(combine_blocks(:,1)==ib)
        bk_raster = [bk_raster raster([raster.block]== combine_blocks(combine_blocks(:,1)==ib,2) )];
        bk_str = [num2str(ib) num2str(combine_blocks(combine_blocks(:,1)==ib,2))];
    else
        bk_str = num2str(ib);
    end
    
    % Convert behavioral state to integer
%     Beh = nan(size(bk_raster));
%     Beh(strcmp({bk_raster.behaving},'P')) = 0;
%     Beh(strcmp({bk_raster.behaving},'A')) = 1;
%     Beh(strcmp({bk_raster.behaving},'D')) = 2;
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');

for ip = 1:max(np)
    
    alldata = bk_raster(np==ip);
    [depths,~,ndpth] = unique([alldata.AMdepth]);

    hF = figure;
    hold on

    % Go through by increasing AMdepth
    for id = 1:max(ndpth)
        subplot(1,max(ndpth),id)
        hold on
        plot_subplot(alldata(ndpth==id),binsize)        
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
    suptitle(alldata(1).stim_str)                     %%%%% change this too
    
    % Get string for savename
    str_pars = strsplit(num2str(LP_HP_dB_rate(ip,:)));
    if numel(depths)==1
        str_dpth = num2str(depths*100);
    else
        str_dpth = 'ROVED';
    end
    
    % Save figure
    an_dir = fullfile(savedir,subject,'^an_plots',session);
    savename = sprintf('%s_%s_FanoFactor_bin%i_ch%i_clu%i_AM%sHz_%sdpth_%sdB_%s-%s_blk%s',...
        subject,session,binsize,channel,clu,str_pars{4},str_dpth,str_pars{3},str_pars{1},str_pars{2},bk_str);

    if ~exist(an_dir,'dir')
        mkdir(an_dir)
    end
    print(hF,'-depsc',fullfile(an_dir,savename))


end
end

% keyboard 

end



function plot_subplot(datablk,binsize)

% Get behaving datapoints
drinking = find(strcmp({datablk.behaving},'D'));
behaving = find(strcmp({datablk.behaving},'A'));
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
FF     = nan(numel(datablk),1);
FFtime = nan(numel(datablk),500);

% Find minimum number of trials for each stimulus
ntc=nan(numel(datablk),1);
for ii = 1:numel(datablk)
    ntc(ii) = max(datablk(ii).y);
end
ent = min(ntc);

for ii = 1:numel(datablk)
    
    data = datablk(ii);
    
    tVec = data.AMonset : binsize : data.stimDur;
    
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

% Set labels
if datablk(1).AMdepth==0
    jitter_labels = repmat({'unmodulated'},1,numel(datablk));
else
    jitter_labels = [datablk.jitter];
end
jitters = str2double(strtok([datablk.jitter],'_'));
nj_max = sum(jitters==mode(jitters));

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
title([num2str(data(1).AMdepth*100) '% depth'])
hAllAxes = findobj(gcf,'type','axes');
if numel(hAllAxes)==1
    ylabel('avg Fano Factor across time during AM')
end


end








