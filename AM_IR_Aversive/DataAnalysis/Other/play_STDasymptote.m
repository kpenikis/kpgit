function play_STDasymptote
% 
% KP, 2019-06-19
% 


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

rng('shuffle')
Iterations = 200;
SampleSize = [5 10 15 20 25 30 40 50 60 80 100];

Unit_std     = nan(1,17);
Unit_std_nrm = nan(1,17);
for iUn = 1:17
    
    istim = 7;
    FRtr = UnitData(iUn).FR_raw_tr(:,istim);
    ntr  = find(isnan([FRtr; nan]),1,'first')-1;
    FRtr = FRtr(1:ntr);
    
    STD_ii_ss = nan(Iterations,length(SampleSize));
    for ss = SampleSize 
        if ss>ntr, break, end
        
        for ii = 1:Iterations
            
            STD_ii_ss(ii,ss==SampleSize) = std(FRtr(randperm(ntr,ss)));
            
        end
        
    end
    
    Unit_std(iUn)     = mean(STD_ii_ss(:,end));
    Unit_std_nrm(iUn) = Unit_std(iUn)/mean(FRtr);
    
    figure;
    plotSpread(STD_ii_ss,'distributionColors','k','showMM',1)
    set(gca,'xtick',1:length(SampleSize),'xticklabels',SampleSize)
    xlabel('Sub-sample size')
    ylabel('Standard deviations')
    title(sprintf('clu %i',UnitData(iUn).Clu))
    
end

% Find max and mins
Unit_std_nrm'

keyboard


end








