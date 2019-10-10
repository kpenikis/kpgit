% call while paused in << plotPopNormMPH >> 
% 

ir=1;

edg_data = zFR_vec(:,[ 1:floor(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):ceil(1000/AMrates(ir)) ],ir);
mid_data = zFR_vec(:, ceil(0.25*1000/AMrates(ir)) : floor(0.75*1000/AMrates(ir)),ir);


zts = [0.1 0.5 0.75 1 1.25 1.5];

for zthresh = zts
    
    zdata(zdata>zthresh) = zthresh;
    zdata(zdata<zthresh) = 0;
    
    
    figure;
    subplot(1,2,1); hold on
    plot(size(edg_data,2)*[-1 1],size(edg_data,2)*[-1 1],'k')
    % plot([-1 1],[0 0],'k')
    % plot([0 0],[-1 1],'k')
    axis square
    
    % edg = nanmean(sign(edg_data) .* double((edg_data)>0.5),2);
    % mid = nanmean(sign(mid_data) .* double((mid_data)>0.5),2);
    edg = sum(edg_data>zthresh,2);
    mid = sum(mid_data>zthresh,2);
    
    plot(edg,mid,'ok')
    plot(edg(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1), mid(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1),'or')
    
    xlim(size(edg_data,2)*[0 1])
    ylim(size(edg_data,2)*[0 1])
    xlabel('half pd around trough')
    ylabel('half pd around peak')
    title(sprintf('# ms above z=%0.2f', zthresh))
    
    
    subplot(1,2,2); hold on
    histogram((edg-mid)./(edg+mid),-1.025:0.05:1.025)
    axis square
    xlim([-1.1 1.1])
    xlabel('(ms trough - ms peak) / sum')
    ylabel('# cells')
    
    zt = char(string(zthresh));
    zt(strfind(zt,'.')) = '-';
    print_eps_kp(gcf,fullfile(savedir,'CatTiming',['Nms_above_' zt]))
    
end


figure;
subplot(1,2,1); hold on
plot([-1 15],[-1 15],'k')
axis square

edg = max(zFR_vec(:,[1:ceil(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):ceil(1000/AMrates(ir))],ir),[],2);
mid = max(zFR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),[],2);

plot(edg,mid,'ok')
plot(edg(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1), mid(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1),'or')

xlim([-1 20])
ylim([-1 20])
xlabel('Edge')
ylabel('Middle')
title('Peak zFR during edges vs middle of MPH')


subplot(1,2,2); hold on
plot([-1 15],[-1 15],'k')
axis square

edg = median(zFR_vec(:,[1:ceil(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):ceil(1000/AMrates(ir))],ir),2);
mid = median(zFR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),2);

plot(edg,mid,'ok')
plot(edg(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1), mid(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1),'or')

xlim([-1 20])
ylim([-1 20])
xlabel('Edge')
ylabel('Middle')
title('Mean zFR during edges vs middle of MPH')



% modeZ = median(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),2);


figure;
edg = max(zFR_vec(:,[1:ceil(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):ceil(1000/AMrates(ir))],ir),[],2);
mid = mode(round(zFR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),1),2);

plot(edg,mid,'ok')
hold on
plot(edg(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1), mid(max(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),[],2) > 1),'or')

xlim([-1 20])
ylim([-1 6])
xlabel('Edges - max zFR')
ylabel('Median of middle')



