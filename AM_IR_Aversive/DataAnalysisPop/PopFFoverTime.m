function PopFFoverTime
% PopFFoverTime
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

%% Load data

close all
fn = set_paths_directories('','',1);

% Load spikes data (created in cumulativeSpikeCount)
load(fullfile(fn.figs,'CmSpCount',['Cell_Time_Trial_Stim_' num2str(20) 'trs']))

%% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');     %[left bottom width height]
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%%
binwin = ones(1,10);

nTrials = size(Cell_Time_Trial_Stim,3);
nStim   = size(Cell_Time_Trial_Stim,4);

for ist = 1:nStim
    
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

