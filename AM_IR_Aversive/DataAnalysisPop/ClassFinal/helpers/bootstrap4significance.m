function [UnSig,ConfIntervals] = bootstrap4significance(CRtable)
% UnSig = bootstrap4significance(CRtable)
%
% 

UnSig = true(size(CRtable,1),1);
Nbs = 10000;

ConfIntervals = nan(size(CRtable,1),2);

for iUn = 1:size(CRtable,1)
    
    Res = CRtable(iUn,:).Results{:};
    
    threshNzeros = round(Nbs*0.05);
    
    dp_bs = nan(Nbs,1);
    for ii = 1:Nbs
        irnd = ceil(rand(size(Res,3),1)*size(Res,3));
        ConfMat = mean(Res(:,:,irnd),3);
        dp_bs(ii)  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);
    end
    
    if sum(dp_bs<=0) >= threshNzeros
        UnSig(iUn) = false;
    end
    
    ConfIntervals(iUn,:) = prctile(dp_bs,[5 95]);
end

% dps = CReach.dprime;
% 
% [newdps,isort] = sort(dps,'ascend');
% iSigs = find(UnSig(isort));
% iNSs = find(~UnSig(isort));
% 
% figure;
% plot(iNSs,newdps(iNSs),'xk')
% hold on
% plot(iSigs,newdps(iSigs),'ok','MarkerFaceColor','k')
% 
% 
% figure;
% plot(find(~UnSig(isort)),dps(isort(~UnSig(isort))),'xk')
% hold on
% plot(find(UnSig(isort)),dps(isort(UnSig(isort))),'ok','MarkerFaceColor','k')
% 

