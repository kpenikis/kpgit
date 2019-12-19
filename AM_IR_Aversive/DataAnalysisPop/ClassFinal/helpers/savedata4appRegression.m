% Load CTTS_AM

% Filter CTTS matrix as in MasterClass (to remove NaNs) 

ETTS=Env_Time_Trial_Stim(theseCells,AnWin,:,theseStim);

TTC_8=[];
for itr = 1:size(CTTS,3)
    TTC_8 = [TTC_8 CTTS(:,:,itr,8)];
end

TTE_8=[];
for itr = 1:size(ETTS,3)
    TTE_8 = [TTE_8 mean(ETTS(:,:,itr,8),1,'omitnan')];
end

InputData = [TTC_8' TTE_8'];
save(fullfile(savedir,'SVRtestdata_8'),'InputData','-v7.3')


% Next, try training on the mean of trials 
