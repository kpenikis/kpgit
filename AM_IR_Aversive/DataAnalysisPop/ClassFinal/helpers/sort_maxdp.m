function [maxdps,iSUdps] = sort_maxdp(CR)


nStim = size(CR.Results{1},1);

% For each cell, get its PC & d' for each stimulus
pcStim = nan(size(CR,1),nStim);
dpStim = nan(size(CR,1),nStim);
% Sparss = nan(size(CR,1),1);

for ii = 1:size(CR,1)
        
    ConfMat = mean(CR(ii,:).Results{:},3,'omitnan');
    
    pcStim(ii,:) = diag(ConfMat)';
    
    ConfMat(ConfMat==0) = 0.01;
    ConfMat(ConfMat==1) = 0.99;
    for ist = 1:size(ConfMat,1)
        othSt = 1:size(ConfMat,1);
        dpStim(ii,ist) =  norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(othSt~=ist,ist)),0,1);
    end
%     Sparss(ii) = calculateSparseness(dpStim(ii,:)');
end


[maxdps,iSUdps] = sort(max(dpStim,[],2),'descend');


end