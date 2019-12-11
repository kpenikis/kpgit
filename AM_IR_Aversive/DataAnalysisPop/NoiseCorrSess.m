function NoiseCorrSess(SUBJECT,SESSION)
% 
% For each pair, plot Signal Corr vs Total Corr
%  Noise corr = total - signal 
%  
%  Signal correlation: combo of tuning and temporal pattern: concatenate PSTHs or MPHs
%    [Cell 1 reshape(stim X time), avg trials] vs [Cell 2 reshape(stim X time), avg trials]
%
%  Total correlation: full matrix, preserving trials; reshape and concatenate
%    [Cell 1 reshape(trials X stim X time)] vs [Cell 2 reshape(trials X stim X time)]
%
% KP, 2019-10
% 

close all

%% Load data

fn = set_paths_directories('','',1);

% Get unit and raster data
[UnitInfo, UnitData, Info, ~, ~, ~, artifactTrs, iUnOrig ] = collectRasterDataSession(SUBJECT,SESSION);

% Load rasters
load(fullfile(fn.figs,'/StimClass/Cell_Time_Trial_Stim_simtrs.mat'))


%%

% Set up gaussian window for smoothing
bs_gaus = 1000;
window = gausswin(bs_gaus);
window = window-min(window);
window = window/sum(window);

% Filter to just this session's cells
Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(iUnOrig,:,:,:);

%
Dur = size(Cell_Time_Trial_Stim,2);

% Get stimuli presented
theseStim = find(permute(~isnan(sum(Cell_Time_Trial_Stim(1,:,1,:),2)),[1 4 2 3]));

% Get minimum number of trials
Ntrs   = nan(numel(theseStim),1);
Ntrs11 = nan(numel(theseStim),1);
for iStim = theseStim
    Ntrs(iStim)   = find(~isnan(Cell_Time_Trial_Stim( 1,1,:,iStim)),1,'last');
    Ntrs11(iStim) = find(~isnan(Cell_Time_Trial_Stim(11,1,:,iStim)),1,'last');
end

%% Calculate

SigCorrs = [];
TotCorrs = [];

for iC1 = 1:(size(Cell_Time_Trial_Stim,1)-1)
    for iC2 = (iC1+1):size(Cell_Time_Trial_Stim,1)
        
        % Gather data for this pair
        SignalC1 = [];
        SignalC2 = [];
        TotalC1  = [];
        TotalC2  = [];
        
        for iStim = theseStim
            
            SignalC1 = [SignalC1 mean(Cell_Time_Trial_Stim(iC1,:,1:min(Ntrs),iStim),3)];
            SignalC2 = [SignalC2 mean(Cell_Time_Trial_Stim(iC2,:,1:min(Ntrs),iStim),3)];
            
            TotalC1 = [TotalC1 reshape(Cell_Time_Trial_Stim(iC1,:,1:min(Ntrs),iStim),1,Dur*min(Ntrs))];
            TotalC2 = [TotalC2 reshape(Cell_Time_Trial_Stim(iC2,:,1:min(Ntrs),iStim),1,Dur*min(Ntrs))];
            
        end
        
        % Convolve data signals
        TotalC1  = conv(TotalC1, window,'same');
        TotalC2  = conv(TotalC2, window,'same');
        
        % Extra step for convolving signal, so no dip at edges
        Sig1 = conv([SignalC1 SignalC1 SignalC1],window,'same');
        Sig2 = conv([SignalC2 SignalC2 SignalC2],window,'same');
        SignalC1 = Sig1(numel(SignalC1)+(1:numel(SignalC1)));
        SignalC2 = Sig2(numel(SignalC2)+(1:numel(SignalC2)));
        
        
        % Calculate signal correlation
        corrSig = corrcoef(SignalC1,SignalC2);
        SigCorrs = [SigCorrs corrSig(1,2)];
        
        % Calculate total correlation
        corrTot = corrcoef(TotalC1,TotalC2);
        TotCorrs = [TotCorrs corrTot(1,2)];
        
    end
end

figure;
plot([-1 1],[-1 1],'k')
hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')

plot(SigCorrs,TotCorrs,'ok')
% xlim([-0.1 0.2])
% ylim([-0.05 0.05])
axis square
xlabel('Signal Correlation')
ylabel('Total Correlation')
title(sprintf('Correlation Relationships, %ims gauss win',bs_gaus))


print_eps_kp(gcf,fullfile(fn.figs,'NoiseCorr',sprintf('CorrRelationships_%ims',bs_gaus)))



% Plot data being correlated, for last pair

figure; hold on
plot(SignalC1,'k')
plot(SignalC2)
xlim([0 numel(SignalC1)])
title('Signal data, last cell pair')
print_eps_kp(gcf,fullfile(fn.figs,'NoiseCorr',sprintf('SignalData_%ims',bs_gaus)))

figure; hold on
plot(TotalC1,'k')
plot(TotalC2)
xlim([0 numel(TotalC1)])
title('Total data, last cell pair')
print_eps_kp(gcf,fullfile(fn.figs,'NoiseCorr',sprintf('TotalData_%ims_zoom',bs_gaus)))



return


% Test reshaping
TestMat = [1:10; 11:20; 21:30];
TestMat(:,:,2) = TestMat(:,:,1)+30;
TestMat(:,:,3) = TestMat(:,:,1)+60;

reshape(TestMat(1,:,:),1,size(TestMat,2)*size(TestMat,3))



end