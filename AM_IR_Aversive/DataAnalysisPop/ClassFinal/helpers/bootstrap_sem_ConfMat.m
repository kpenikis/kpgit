function sem = bootstrap_sem_ConfMat(ConfMat,yaxvalue,correction)

Nbs = 500;
Sbs = size(ConfMat,3);
rng('shuffle')

ConfMat(ConfMat==0) = correction;
ConfMat(ConfMat==1) = 1-correction;

switch yaxvalue
    
    case 'dp'
        dps = nan(Nbs,size(ConfMat,1));
        for ibs = 1:Nbs
            thisCM = mean(ConfMat(:,:,randi(Sbs,[1 1])),3);
        for ist = 1:size(ConfMat,1)
            dps(ibs,ist) = norminv(thisCM(ist,ist),0,1) - norminv(mean(thisCM(setdiff(1:end,ist),ist)),0,1);
        end
        end
        sem = std(dps,1)./sqrt(Nbs);
        
    case 'PC'
        pcs = nan(Nbs,size(ConfMat,1));
        for ibs = 1:Nbs
            pcs(ibs,:) = diag(mean(ConfMat(:,:,randi(Sbs,[1 1])),3))';
        end
        sem = std(pcs,1)./sqrt(Nbs);        
end

end