

plotdata  = ThisData(i_sorted,:,ist).*1000;
xtickset  = [0 size(plotdata,2)];
xticklabs = xtickset;
ndp = size(plotdata,1);


% Render plot
imagesc(log10(plotdata(1:ndp,:)))
%         caxis([0 log10(max(Boundaries))+0.25])
caxis([0 2.5])
cmocean('-gray')


% Finish plot
xlim([1 size(plotdata,2)])
ylim([0.5 ndp+0.5])
set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
set(gca,'xtick',xtickset,'xticklabel',xticklabs)

if exist('yticklab','var')
    set(gca,'ytick',ytickset,'yticklabel',yticklab)
else
    set(gca,'ytick',[])
end
if ist==1
    ylabel(sortBy)
end

title(['stim' num2str(ist)])
box off
axis fill