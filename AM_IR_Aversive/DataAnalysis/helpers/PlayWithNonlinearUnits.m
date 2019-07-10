
fn  = set_paths_directories;
nUn = 5;

% Load Unit files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load PredObs FR data
q = load(fullfile(fn.figs,'PredObs','Pdata_FR.mat'));
POdata = q.Pdata;
clear q

% Load PredObs TRF data
q = load(fullfile(fn.figs,'TRF','Pdata_TRF.mat'));
TRFdata = q.Results;
clear q

% Load Classifier results
q = load(fullfile(fn.processed,'MPHclassifier','ClassData.mat'));
ClassData = q.Data;
clear q


%% PredObs FR off of unity

[POd,iPOd] = sort(abs([POdata.obsIR] - [POdata.predIR]),'descend');
POds = sign([POdata.obsIR] - [POdata.predIR]);
POd = POd.*POds(iPOd);

POdata(iPOd(1:nUn+6),:)

MPH_plot(124) % MPHs no big diffs; see rasters for expl.
              % DB lower than predicted; possibly units fades toward end of recording
              % AC higher than predicted; same reason: later Pdc trials pulling down prediction
              % sooo, if only use pdc trials right around that IR, should be
              % no difference
MPH_plot(126) % Same explanation; FR throughout session severely decays
              % TRUE FOR ALL CH 4 THIS SESSION

MPH_plot(73)  % AC lower than predicted; looks like 4 and 8 Hz are suppressed following faster rates

MPH_plot(84)  % DB lower than predicted; not clear why
ClassData(84,2).Res_L1o.dprime 

MPH_plot(123) % see above

MPH_plot(143) % DB lower than predicted; 32Hz Pdc is facilitating

MPH_plot(43)  % n.s. 

MPH_plot(123) % see above

MPH_plot(77)  % DB lower than predicted; troughs of 4 not so low in Pdc; same unit as 73



%% TRF Pred off of unity

[TRFd,iTRFd] = sort(abs([TRFdata.C100_Irr] - [TRFdata.C100_Pdc]),'descend');
TRFds = sign([TRFdata.C100_Irr] - [TRFdata.C100_Pdc]);
TRFd = TRFd.*TRFds(iTRFd);


TRFdata(iTRFd(1:nUn+6),:)
% NO OVERLAP IN UNITS

MPH_plot(121) % weird onset/offset responses

MPH_plot(56)  % 2hz Pdc missing onset spike

MPH_plot(57)  % Pdc: 2hz has onset spike; IR: 32hz spike

MPH_plot(23)  % pretty unreliable

MPH_plot(58)  % sparse/unreliable

MPH_plot(64)  % 



%% Finding highest dprimes

for iUn = 1:size(ClassData,1)
    try
    if any(ClassData(iUn,1).Res_L1o.dprime(:,2) > 1.5)
       disp(iUn) 
    end
    end
end

% Mar30-AM sound drops out! After removing these trials will dprime still be high??
MPH_plot(149)                 % good responses; many high dprimes, but MPHs don't always look v different
TRFdata([TRFdata.iUn]==149,:) % linear
POdata([POdata.iUn]==149,:)   % linear

MPH_plot(118)                 % Very sparse
TRFdata([TRFdata.iUn]==118,:) % IR is a bit better
POdata([POdata.iUn]==118,:)   % linear

MPH_plot(69)                  % good responses; 32hz real differences, 2hz due to previous pd response lingering
TRFdata([TRFdata.iUn]==69,:)  % linear
POdata([POdata.iUn]==69,:)    % linear



%% Correlation of ranks

POdiffUnique = nan(numel(UnitData),2);
for ii = 1:size(POdata,1)
     POdiffUnique(POdata.iUn(ii),POdata.stid(ii)-6) = POdata.obsIR(ii) - POdata.predIR(ii);
end
[foo,iPOdU] = sort(abs(mean(POdiffUnique,2,'omitnan')),'descend');

UnitsN = 1:numel(UnitData);

figure; 
plot((iPOdU),(iTRFd),'.k')

figure;
plot((iPOdU),iTRFd(iPOdU),'.k')




