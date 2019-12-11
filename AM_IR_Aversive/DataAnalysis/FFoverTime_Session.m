function FFoverTime_Session(SUBJECT,SESSION)
% PopFFoverTime
%  
%  Compare variability of avg activity across units (zhou)
%  to variability of each neuron across trials
%
%  Looking at this because of Zhou Wang 2012. They showed FF encoding delta
%  amplitude, when FF calculated from avg PSTHs across all nonmonotonic
%  cells (n=95?). 
%  Here I look at same thing but during individual trials, across ~255
%  cells. Some structure at amplitude transitions begins to pop out in the
%  median FF across trials, but not in individual ones. So to see this
%  degree of structure in single trials, you'd need to be decoding from 
%  ~20x255, or 5100 cells simultaneously.
%  Dependence on binning window: 20 ms too long; 5 ms too short.
%
%  KP, 2019-10
%

minTrs = 15;
binwin = ones(1,50);

%% Load data

% close all
fn = set_paths_directories('','',1);

% Load spikes data (created in cumulativeSpikeCount)
load(fullfile(fn.figs,'StimClass',['Cell_Time_Trial_Stim_simtrs']))

% Load unit files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Filter to just this session
iUnOrig  = find(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitData = UnitData(iUnOrig);
UnitInfo = UnitInfo(iUnOrig,:);

Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(iUnOrig,:,:,:);

theseStim = find(~isnan(Cell_Time_Trial_Stim(1,1,1,:)))';
nTrs = sum(~isnan(sum(permute(Cell_Time_Trial_Stim(1,1,:,theseStim),[3 4 1 2]),2)));
if nTrs<minTrs
    keyboard
end

Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(:,:,1:nTrs,theseStim);

% Convolve each time series
tic
for iUn = 1:size(Cell_Time_Trial_Stim,1)
    for iTr = 1:size(Cell_Time_Trial_Stim,3)
        for ist = 1:numel(theseStim)
            Cell_Time_Trial_Stim(iUn,:,iTr,ist) = conv(Cell_Time_Trial_Stim(iUn,:,iTr,ist),binwin,'same');
        end
    end
end
toc


%% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%%

for ist = 7:numel(theseStim)
    
    FF_popmean = var(mean(Cell_Time_Trial_Stim(:,:,:,ist),3),[],1) ./ mean(mean(Cell_Time_Trial_Stim(:,:,:,ist),3),1);
    
    
    FF_Trs_cell = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,2));
    r_to_PopFF  = nan(size(Cell_Time_Trial_Stim,1),2);
    
    for iUn = 1:size(Cell_Time_Trial_Stim,1)
        
       FF_Trs_cell(iUn,:) = var(Cell_Time_Trial_Stim(iUn,:,:,ist),[],3)  ./ mean(Cell_Time_Trial_Stim(iUn,:,:,ist),3);
       
%        if sum(isnan(FF_Trs_cell(iUn,:)))<250
           rfoo = []; pfoo = [];
           [rfoo,pfoo] = corrcoef(FF_Trs_cell(iUn,:), FF_popmean);
           r_to_PopFF(iUn,:) = [rfoo(1,2) pfoo(1,2)];
%        end
       
    end
    
    figure; hold on
    plot(FF_Trs_cell')
    plot(FF_popmean,'k','LineWidth',6)
    plot(mean(FF_Trs_cell,1,'omitnan'),'y','LineWidth',3)
    
    keyboard
    
    FF = nan(nTrials,size(Cell_Time_Trial_Stim,2));
    
    for it = 1:nTrials
        
        binned_data = [];
        for iUn = 1:size(Cell_Time_Trial_Stim,1)
            binned_data = [binned_data; conv(Cell_Time_Trial_Stim(iUn,:,it,ist),binwin)];
        end
        binned_data = binned_data(:,1:size(Cell_Time_Trial_Stim,2));
        
        FF(it,:) = var(binned_data,1)./mean(binned_data,1);
        
    end
    
    FFmed = median(FF,1);
    
    figure; 
    set(gcf,'Position',widescreen)
    plot(FF')
    hold on
    plot(FFmed,'k','LineWidth',4)
    title(num2str(ist))
    ylim([0.8 1.4])
    
end %ist



end

