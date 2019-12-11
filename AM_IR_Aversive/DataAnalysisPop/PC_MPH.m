function PC_MPH
% 
% PCA population response trajectories for each MPH.
%
% KP, 2019-12
% 

close all

% Options
plotWhich  = 'sp';
theseStim  = 2:4;
bs_gaus    = 5;
spshft     = ''; %_noshift


%% Prepare

fn = set_paths_directories('','',1);

% Load spiking data for population
load(fullfile(fn.figs,'StimClass/Cell_Time_Trial_Stim_simtrs.mat'))
% more stimulus excerpts than just periodic rates; 2:6 correspond to AM rates

% Load pre-made MPH data (smoothed + averaged activity for each cell, for
% one period of each AM rate)
load(fullfile(fn.figs,'PopMPH',['MPHdata' spshft '.mat']))


% Fig settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

colors = [  84  24  69;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 230 120]./255;


%% Gather data

% Set up gaussian window for smoothing
window = gausswin(bs_gaus);
window = window-min(window);
window = window/sum(window);

% Info about data mat
Ncells = size(Cell_Time_Trial_Stim,1);
Dur    = size(Cell_Time_Trial_Stim,2);

% Set stimuli 
nStim = numel(theseStim);
AMrates   = [2 4 8 16 32];


% Calculate average PSTH for each cell/stim
Cell_Time_Stim = nan(Ncells,Dur,nStim);
for iUn = 1:Ncells
    for ist = 1:nStim
        
        thisDur = floor(1000./AMrates(ist));
        
        Cell_Time_Stim(iUn,1:thisDur,ist) = conv(mean(Cell_Time_Trial_Stim(iUn,1:thisDur,:,theseStim(ist)),3,'omitnan'),window,'same');
        
    end
end


% Remove cells with not enough trials
rmv_cells = [];
for iUn = 1:size(FR_vec,1)
    if all(isnan(FR_vec(iUn,:,1)))
        rmv_cells = [rmv_cells iUn];
    end
end

FR_vec(rmv_cells,:,:)          = [];
zFR_vec(rmv_cells,:,:)         = [];
Cell_Time_Stim(rmv_cells,:,:)  = [];


%% PCA and plot result

ii = [1 cumsum(floor(1000./AMrates(1:nStim)))];

% PCA of data from all stimuli

Data     = [];
for ist = 1:nStim % skip 32 hz
    
    thisDur = floor(1000./AMrates(ist));
    
    % Choose to use either spiking data, PSTHs, or z-scored activity
    switch plotWhich
        case 'FR'
            Data = [Data FR_vec(:,1:thisDur,ist)];
        case 'zFR'
            Data = [Data zFR_vec(:,1:thisDur,ist)];
        case 'sp'
            Data = [Data Cell_Time_Stim(:,1:thisDur,ist)];
    end
end

[coefs,~,latent] = pca(Data);

fprintf('\nPC1 %0.1f%%, PC2 %0.1f%%, PC3 %0.1f%%\n\n',...
    100*latent(1)/sum(latent), 100*latent(2)/sum(latent), 100*latent(3)/sum(latent) )


% Plot 

hf=figure; hold on

for ist = 1:nStim % skip 32 hz
    
    pkEidx = floor((ii(ist+1)-ii(ist))/2) + ii(ist);
    pkRidx = floor((ii(ist+1)-ii(ist))/4) + ii(ist);
    
    ip(ist)=plot3( coefs(1+ii(ist):ii(ist+1),1), coefs(1+ii(ist):ii(ist+1),2), coefs(1+ii(ist):ii(ist+1),3), '-','Color',colors(ist,:),'LineWidth',0.5);
    ip(ist)=plot3( coefs(1+ii(ist):ii(ist+1)), coefs(1+ii(ist):ii(ist+1),2), coefs(1+ii(ist):ii(ist+1),3), '.','Color',colors(ist,:),'MarkerSize',10);
    plot3( coefs(1+ii(ist),1), coefs(1+ii(ist),2), coefs(1+ii(ist),3), 'o','MarkerFaceColor',colors(ist,:),'MarkerSize',20,'MarkerEdgeColor','k')
    plot3( coefs(pkEidx,1), coefs(pkEidx,2), coefs(pkEidx,3), 'o','MarkerEdgeColor',colors(ist,:),'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',2)
%     plot3(coefs(pkRidx,1),coefs(pkRidx,2),coefs(pkRidx,3),'o','MarkerEdgeColor',colors(ist,:),'MarkerFaceColor','c','MarkerSize',10,'LineWidth',2)
    
end

axis square
xlabel(sprintf('PC1, %0.1f%%',100*latent(1)/sum(latent)))
ylabel(sprintf('PC2, %0.1f%%',100*latent(2)/sum(latent)))
zlabel(sprintf('PC3, %0.1f%%',100*latent(3)/sum(latent)))

legend(ip,cellfun(@(x) sprintf('%i Hz',x),num2cell(AMrates(1:nStim)),'UniformOutput',false),'Location','best')

title(sprintf('AM periods in PC space, %s',plotWhich))





print_eps_kp(hf,fullfile(fn.figs,'PCtraj',sprintf('PC_MPH_%s%s',plotWhich,spshft)))



end



