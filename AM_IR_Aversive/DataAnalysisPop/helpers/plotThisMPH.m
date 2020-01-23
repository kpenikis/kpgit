% 
data = ThisData(theseCells,1:ceil(1000/AMrates(ir)),ir);

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Sort cells
if ~exist('i_sorted','var')
    [i_sorted,sortdata] = sort_thLat(data);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Get ready to plot

plotdata  = data(i_sorted,:);
xtickset  = [0 ceil(1000/AMrates(ir))/2 ceil(1000/AMrates(ir))];
xticklabs = {'0' 'pi' '2*pi'};
ndp = size(plotdata,1); %sum(sum(isnan(plotdata),2)==0);

if clipZ>0
    plotdata(plotdata>clipZ) = clipZ;
    plotdata(plotdata<clipZ) = 0;
end

% Render plot
switch useFR
    case 'z'
        imagesc(plotdata(~any(isnan(plotdata),2),:))
        caxis([min(Boundaries) max(Boundaries)])
        cmocean('balance','pivot',0) %curl
    case 'log'
        imagesc(log10(plotdata(1:ndp,:)))
%         caxis([0 log10(max(Boundaries))+0.25])
        caxis([0 1.75])
        cmocean('-gray')
end


if labelNS
    % Add markers to label NS cells
    hold on
    plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
    plot(ceil(1000/AMrates(ir)),find(flagNS),'.','Color',[0.01 0.57 0.44])
end


% Finish plot
xlim([0 ceil(1000/AMrates(ir))])
ylim([0.5 ndp+0.5])
set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
set(gca,'xtick',xtickset,'xticklabel',xticklabs)

if ~sortHere
    BoundMat = cumsum(sortdata(i_sorted(1:ndp),1)>=Boundaries);
    ytplc = []; ytlab = [];
    for ib = numel(Boundaries):-1:1
        yUn = find(BoundMat(:,ib)==1,1,'first');
        if ~ismember(yUn,ytplc)
            ytplc = [ytplc yUn];
            ytlab = [ytlab Boundaries(ib)];
        end
    end
    set(gca,'ytick',fliplr(ytplc),'yticklabel',fliplr(ytlab))

else 
    if exist('yticklab','var')
        set(gca,'ytick',ytickset,'yticklabel',yticklab)
    else
        set(gca,'ytick',[])
    end
end

% CluLabels = [UnitData(theseUnits(i_sorted)).Clu];
% set(gca,'ytick',theseUnits,'yticklabel',CluLabels)

title([num2str(AMrates(ir)) ' Hz'])
box off
axis fill