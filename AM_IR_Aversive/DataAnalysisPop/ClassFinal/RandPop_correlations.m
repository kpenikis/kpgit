

% Load CR random pools
% and associated CTTS


% For each group of N (10) cells, for each stimulus,
% estimate projection measure and compare to d' (just that stim)
nstim = 8;
ncell = size(CTTS_AM,1);
ncPop = mode(CR.nC);

nbs   = 100;
rng('shuffle')

ProjEsts = nan(nstim,size(CR,1),ncPop);
dpStims  = nan(nstim,size(CR,1));

nTrialMat = nan(size(CTTS_AM,1),size(CTTS_AM,4));
for ist = 1:size(CTTS_AM,4)
    CT  = permute(sum(CTTS_AM(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

for ipop = 1:size(CR,1)
    
    iCs  = CR(ipop,:).SUids{:};
    Res  = mean(CR(ipop,:).Results{:},3); 
    
    
    for istim = 1:nstim
        
        ntr = nTrialMat(iCs,istim);
        
        theseProjs = nan(numel(iCs),nbs);
        for i = 1:nbs
            
            theseTrs = ceil(rand(numel(iCs),11).*ntr);
            
            for ic = 1:numel(iCs)
                theseProjs(ic,i) = corr(mean(CTTS_AM(iCs(ic),:,theseTrs(ic,1:10),istim),3)', CTTS_AM(iCs(ic),:,theseTrs(ic,11),istim)');
            end
            
        end
        
        intProjEst =  mean(theseProjs,2,'omitnan');
        intProjEst(sum(isnan(theseProjs),2)>(0.5*nbs)) = nan;
        
        ProjEsts(istim,ipop,:) = intProjEst;
        
    end %istim
    
    dpStims(:,ipop) = dp_from_ConfMat(Res,0.01);
    
    
end %ipop


figure;
for istim = 1:nstim
    
    subplot(2,4,istim)
    hold on
%     plot(max(ProjEsts(istim,:,:),[],3),dpStims(istim,:),'.k')
    plot(median(ProjEsts(istim,:,:),3,'omitnan'),dpStims(istim,:),'.k')
    
    xlim([0 0.2])
    ylim([0 4])
end


figure;
plot(median(max(ProjEsts,[],3),1),mean(dpStims,1),'.k')







