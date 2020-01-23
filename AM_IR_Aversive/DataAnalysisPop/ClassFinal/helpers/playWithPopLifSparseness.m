
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<0.43);% & [UnitData(theseCells).BaseFR]'>2);

[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');
iBest5  = iRS(iSUdps(1:5));
iBest10 = iRS(iSUdps(1:10));
iBest15 = iRS(iSUdps(1:15));

[dps,iSUdps] = sort(CReach(iNS,:).dprime,'descend');
iBest5  = iNS(iSUdps(1:5));

newcolors = cmocean('thermal',5);


PSTHs     = mean(CTTS,3,'omitnan');

PopSps = zeros(size(PSTHs,4),size(PSTHs,2));

figure; 
hold on
for ist = 1:size(PSTHs,4)
    for ims = 3:(size(PSTHs,2)-2)
        PopSps(ist,ims) = calculateSparseness(sum(PSTHs(:,ims+(-2:2),:,ist),2));
    end
   
    subplot(4,2,ist);
    hold on
    yyaxis left
%     plot(PopSps(ist,:),'Color',0.75*[1 1 1],'LineWidth',8)
    fill([1:size(PopSps,2) size(PopSps,2) 1 1],...
        [PopSps(ist,:) 0 0 PopSps(ist,1)],...
        0.8*[1 1 1], 'EdgeColor',0.5*[1 1 1])
    xlim([0 size(PopSps,2)])
    ylim([0.5 0.8])
    
    % Best
    yyaxis right
    set(gca,'ColorOrder',newcolors);
    
    % each best unit
    plot(PSTHs(iBest5,:,:,ist)'.*1000,'-','LineWidth',2)    
    ylim([0 600])
    
    % sum best & sum rest
%     plot(sum(PSTHs(iRS(iSUdps(21+(0:30))),:,:,ist),1),'-r','LineWidth',2)
%     plot(sum(PSTHs(iBest10,:,:,ist),1),'-b','LineWidth',2)
    
    % ratios
%     plot(mean(PSTHs(iRS(iSUdps(11:15)),:,:,ist),1)./mean(PSTHs(iRS(iSUdps(16:end)),:,:,ist),1),'-b','LineWidth',2)
%     plot(mean(PSTHs(iRS(iSUdps(6:10)),:,:,ist),1)./mean(PSTHs(iRS(iSUdps(16:end)),:,:,ist),1),'-g','LineWidth',2)
%     plot(mean(PSTHs(iBest5,:,:,ist),1)./mean(PSTHs(iRS(iSUdps(16:end)),:,:,ist),1),'-m','LineWidth',2)
end


print_eps_kp(gcf,fullfile(rootdir,'Best5_NS_PopSps'))





PSTHs     = mean(CTTS,3,'omitnan');
NormPSTHs = PSTHs ./ sum(sum(PSTHs,2),4);

% Lifetime sparseness 
LifSps = nan(size(PSTHs,1),1);

for iun = 1:size(PSTHs,1)
    
    CatResp = reshape(permute(PSTHs(iun,:,:,:),[2 4 1 3]),[size(PSTHs,2)*size(PSTHs,4) 1]);
    for ims = 1:size(PSTHs,2)
        LifSps(iun,1) = calculateSparseness(CatResp);
    end
end

WeightedResps = nan(size(PSTHs,4),size(PSTHs,2));
for ist = 1:size(PSTHs,4)
    
    WeightedResps(ist,:) = mean(LifSps.*PSTHs(:,:,:,ist),1);
    
end
WeightedResps = WeightedResps/max(max(WeightedResps)); 

WtNormResps = nan(size(NormPSTHs,4),size(NormPSTHs,2));
for ist = 1:size(NormPSTHs,4)
    
    WtNormResps(ist,:) = mean(LifSps .* NormPSTHs(:,:,:,ist) ,1);
    
end
WtNormResps = WtNormResps/max(max(WtNormResps)); 



% Population sparseness
PopSps = nan(size(PSTHs,4),size(PSTHs,2));

figure; 
hold on
for ist = 1:size(PSTHs,4)
    
    for ims = 1:size(PSTHs,2)
        PopSps(ist,ims) = calculateSparseness(PSTHs(:,ims,:,ist));
    end
   
    subplot(4,2,ist);
    hold on
    plot(PopSps(ist,:),'LineWidth',4)
    xlim([0 size(PopSps,2)])
    ylim([0.2 1])
    
    plot(mean(PSTHs(:,:,:,ist),1)./max(max(mean(PSTHs,1),[],2)),'m','LineWidth',2)
%     
%     % Add weighted lifetime sparseness
%     plot(WeightedResps(ist,:),'LineWidth',2)
%     plot(WtNormResps(ist,:),'LineWidth',2)
    
end







