function ap_plot_FRresp(subject,session,channel,clu,raster)
% Plots the average firing rate during the AM portion of a stimulus. 

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

if nargin<5
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end

% Find stimuli with more than 8 trials
raster = raster(cellfun(@length,{raster.tr_idx}) > 8);

% Set labels
jitters = str2double(strtok([raster.jitter],'_'));
jitter_labels = [raster.jitter];
nj_max = sum(jitters==mode(jitters));

% Get FR and std of FR for error bars
FRvec = [raster.nSpk];
FRstd = nan(1,numel(raster)); FRtr=struct();
for ir = 1:numel(raster)
    x = raster(ir).x;
    y = raster(ir).y;
    FR_t = nan(1,max(y));
    for it = 1:max(y)
        FR_t(it) = numel(x(x(y==it) > raster(ir).AMonset & x(y==it) < raster(ir).stimDur)) / ((raster(ir).stimDur-raster(ir).AMonset)/1000);
    end
    FRtr(ir).tr = FR_t;
    FRstd(ir) = std(FR_t,'omitnan');
end


baselineFR = sum([raster.x]<0) / (raster(1).window_ms(1)/-1000) / sum(cellfun(@(x) ( numel(x) ), {raster.tr_idx}));
hF = figure; hold on

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
    hold on
    errorbar(ii:ii+nj_max-1, ybarvec(ii:ii+nj_max-1), ebarvec(ii:ii+nj_max-1),...
        'LineStyle','none', 'Color',[0 0 0] , 'LineWidth',1 )

end

mean_j0FR = mean(ybarvec(xbarvec==0),'omitnan');

% Figure properties
set(gca,'xtick',1:length(xbarvec),'xticklabel',xbarlab,'TickLabelInterpreter', 'none','XTickLabelRotation',45)
plot([0 length(xbarvec)+1],[mean_j0FR mean_j0FR],'--b')
set(gca,'xlim',[0 length(xbarvec)+1])
title(strtok(raster(1).stim_str,','))
ylabel('avg FR during AM')
hold off

% Save figure
an_dir = fullfile(savedir,subject,'^an_plots',session);
savename = sprintf('%s_%s_FRresp_ch%i_clu%i',subject,session,channel,clu);
if ~exist(an_dir,'dir')
    mkdir(an_dir)
end
print(hF,'-depsc',fullfile(an_dir,savename))


end

