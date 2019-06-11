
% Plot MPH context classification results 
% separately for each putative cell type


thinData = [];
thickData = [];

figure;
hold on
for iUn = 1:numel(UnitData)
    
    plot(UnitInfo(iUn,:).TroughPeak,    Data(iUn,2).Res_L1o.dprime(:,2)',...
        'ok','LineWidth',2)
    
    if UnitInfo(iUn,:).TroughPeak < 0.45
        thinData = [thinData Data(iUn,2).Res_L1o.dprime(:,2)'];
    else
        thickData = [thickData Data(iUn,2).Res_L1o.dprime(:,2)'];
    end
    
end


figure;
histogram(thinData,0:0.1:4,'Normalization','probability')
hold on
histogram(thickData,0:0.1:4,'Normalization','probability')

