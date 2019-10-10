


% Remove empty units
iUn_GoodData = find(sum(~isnan(zFR_vec(:,:,ir)),2)==ceil(1000/AMrates(ir)));
GoodData = zFR_vec(iUn_GoodData,1:ceil(1000/AMrates(ir)),ir);


[COEFF, SCORE, LATENT] = pca(GoodData');

% figure;
% plot(GoodData(:,ceil(end/6)),GoodData(:,ceil(end/2)),'ok')
% hold on
% plot(SCORE(ceil(end/6),:),SCORE(ceil(end/2),:),'ob')


% Calulate linkages
Y = pdist(SCORE');     % cells X time
Z = linkage(Y,'ward');

% leafOrder = optimalleaforder(Z,Y);

% Plot dendrogram
figure;
set(gcf,'Position',fullscreen)

subplot(1,6,2:6);      % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'ColorThreshold',10,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
set(gca,'tickdir','out')

RespCluLabels = cluster(Z,'maxclust',10);

% Plot MPH
subplot(1,6,1);

imagesc( SCORE(:,fliplr(outperm))' )
% imagesc( GoodData(fliplr(outperm),:) )

caxis([-clipVal clipVal])
cmocean('balance','pivot',0)
xlim([0 ceil(1000/AMrates(ir))])
ylim([0 size(GoodData,1)])
set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
box off
axis fill

suptitle([num2str(AMrates(ir)) ' Hz'])


