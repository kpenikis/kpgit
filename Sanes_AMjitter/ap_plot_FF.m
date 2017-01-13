function ap_plot_FF(subject,session,channel,clu,raster)
% Plots a folded histogram of spiking during the AM portion of a stimulus.
% Note that periods are not all of equal length so n trials tapers off.
% CV = std / mean;
% FF = var / mean;

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

if nargin<5
    savedir  = '/Users/kpenikis/Documents/SanesLab/Data/processed_data';
    savename = sprintf('%s_sess-%s_raster_ch%i_clu%i',subject,session,channel,clu);
    load(fullfile(savedir,subject,savename))
end


jitters = str2double(strtok([raster.jitter],'_'));

% Set up for creating histo
binsize = 100;

FF     = nan(numel(raster),1);
FFtime = nan(numel(raster),500);

ntc=nan(numel(raster),1);
for ii = 1:numel(raster)
    ntc(ii) = max(raster(ii).y);
end
ent = min(ntc);

for ii = 1:numel(raster)
    
    data = raster(ii);
    
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

% Now set the order for plotting
jitters = str2double(strtok([raster.jitter],'_'));
nj_max = sum(jitters==mode(jitters));
jitter_labels = [raster.jitter];

hF = figure;

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
title(strtok(raster(1).stim_str,','))
ylabel('avg Fano Factor across time during AM')
hold off

% Save figure
an_dir = fullfile(savedir,subject,'^an_plots',session);
% savename = sprintf('%s_%s_FanoFactor_bin%i_ch%i_clu%i_AM%iHz_jitter%s_%idpth_%idB_%i-%i_%s_blk%i',...
%     subject,session,binsize,channel,clu,...
%     data.AMrate, data.jitter, round(data.AMdepth*100),...
%     data.dB, data.HP, data.LP, data.behaving, data.block);
savename = sprintf('%s_%s_FanoFactor_bin%i_ch%i_clu%i',subject,session,binsize,channel,clu);
if ~exist(an_dir,'dir')
    mkdir(an_dir)
end
print(hF,'-depsc',fullfile(an_dir,savename))



end







