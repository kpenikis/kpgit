function exploreNoiseCorr(SUBJECT,SESSION)
% 
% Must be single session.
% Based on playNoiseCorr.
% After plotting PopMPHs and talk 1 with Cristina.
% 
% KP, 2019-10
% 


global AMrates useFR Boundaries trMin
close all
%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
%~~~~~~~~~~~~~~~~~~~~~

%% Load data
fn = set_paths_directories('','',1);

% Get unit and raster data
[UnI, UnD, Info, TrialData, Clusters, StimResp, ~, iUnOrig ] = collectRasterDataSession(SUBJECT,SESSION);

% filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


% Load PopMPH data
savedir = fullfile(fn.figs,'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
load(fullfile(savedir,'MPHdata.mat'))

FR_vec   = FR_vec(iUnOrig,:,:); 
FR_Warn  = FR_Warn(iUnOrig,:,:); 
zFR_vec  = zFR_vec(iUnOrig,:,:); 
zFR_Warn = zFR_Warn(iUnOrig,:,:); 

if size(FR_vec,1) ~= size(UnI,1)
    warning('saved data files have different n units')
    keyboard
%     plotPopNormMPH(1)
%     load(fullfile(savedir,'MPHdata.mat'))
end


% Data settings

switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5 1.75]) ];
end

        
%%

stid = 3;
data = ThisData(:,1:ceil(1000/AMrates(stid-1)),stid-1);

% [i_sorted,sortdata] = sort_thLat(data);
i_sorted = 1:size(data,1);

figure; 
plot(repmat(1:size(data,1),size(data,2),1)*20 + data(flipud(i_sorted),:)','k')


% Choose a cell
iUn_1 = iUnOrig(i_sorted(end-2));
UnitInfo(iUn_1,:);                    %clu 721


% Choose 2 similar
iUn_2 = iUnOrig(i_sorted(end-7));
UnitInfo(iUn_2,:);                    %clu 424

iUn_3 = iUnOrig(i_sorted(end-9));
UnitInfo(iUn_3,:);                    %clu 760

% Choose 2 different
iUn_4 = iUnOrig(i_sorted(end));
UnitInfo(iUn_4,:);                    %clu 879

iUn_5 = iUnOrig(i_sorted(end-16));
UnitInfo(iUn_5,:);                    %clu 699


%% Get noise and signal correlations

% Signal + each Noise, for each pair
% Signal vs avg(Noise) for all pairs

Clu1 = 721;
Clu2 = 879;

Un1  = [UnI.Clu]==Clu1;
Un2  = [UnI.Clu]==Clu2;
spiketimes1 = unique(round(Clusters([Clusters.clusterID]==Clu1).spikeTimes*1000));
spiketimes2 = unique(round(Clusters([Clusters.clusterID]==Clu2).spikeTimes*1000));


%% Now collect R noise values during each stimulus

% Get all stimuli presented with these parameters, given a
% sufficient number of trials without diruptive artifact
% while the animal was drinking
[all_TDidx,Ntrials] = get_clean_trials( TrialData,...
    unique([Info.artifact(UnD(Un1).Channel).trials; Info.artifact(UnD(Un2).Channel).trials]),...
    UnD(Un1).spl,UnD(Un1).lpn,1);
allStim = unique(TrialData.trID(all_TDidx))';

minTrs = min(Ntrials(~isnan(Ntrials)));

ymaxval = 0;

Rnoise_stim       = nan(1,numel(allStim));

FR1_stim          = nan(1,numel(allStim));
FR2_stim          = nan(1,numel(allStim));


for ist = 1:numel(allStim)
    
    stid = allStim(ist);
    
    %% Collect trial indices and timestamps
    
    if stid==3 || stid==6
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
        % Find Pdc trials that follow same rate during ITI
        TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
        
        TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
        TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        
    else
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
        TDidx_iti = [];
    end
    
    
    % Get timestamps of onsets and offsets
    clear t2 t3 Duration t_win
    t2 = TrialData.onset(TDidx);
    t3 = TrialData.offset(TDidx);
    
%     Duration = mode(diff([t2 t3],1,2));
    Duration = 1000;
    
    % Add ITI trials (shortened to match duration)
    if ~isempty(TDidx_iti)
        t2 = [t2; TrialData.onset(TDidx_iti)];
        TDidx = [TDidx; TDidx_iti];
    end
    
    kt     = randperm(length(t2),min([minTrs trMin*2]));
    t2     = t2(kt);
    TDidx  = TDidx(kt);
    
    t3     = t2 + Duration;
    
    
    %%
    FR1 = [];
    FR2 = [];
    
    for it = 1:numel(TDidx)
        
        FR1 = [ FR1 sum( spiketimes1>t2(it) & spiketimes1<=t3(it) ) ];
        FR2 = [ FR2 sum( spiketimes2>t2(it) & spiketimes2<=t3(it) ) ];
        
    end %it
    
    r = corrcoef(FR1,FR2);
    Rnoise_stim(ist) = r(1,2);
    
    FR1_stim(ist) = mean(FR1);
    FR2_stim(ist) = mean(FR2);
    
end %ist

% Signal correlation
r = corrcoef(FR1_stim,FR2_stim);
Rsig = r(1,2);


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

% Set colors
colors = [ 150 150 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];
        

%%

hf=figure; hold on
set(hf,'Position',widescreen)

for ist = 1:numel(allStim)
    
    subplot(2,numel(allStim),ist)
    if ist>1 && ist<7
        plot(ThisData(i_sorted(Un1),1:ceil(1000/AMrates(ist-1)),ist-1),'k','LineWidth',3)
        hold on
        plot(ThisData(i_sorted(Un2),1:ceil(1000/AMrates(ist-1)),ist-1),'b','LineWidth',3)
        xlim([0 ceil(1000/AMrates(ist-1))])
        ylim([0 50])
    elseif ist==1
        plot(FR_Warn(i_sorted(Un1),:),'k','LineWidth',3)
        hold on
        plot(FR_Warn(i_sorted(Un2),:),'b','LineWidth',3)
        xlim([0 500])
        ylim([0 50])
        ylabel('FR (sp/s)')
        xlabel('ms')
    end
    set(gca,'Color','none')
    
    % Correlations
    subplot(2,numel(allStim),ist+numel(allStim))
    plot([0 1],[0 0],'--k','LineWidth',2)
    hold on
    
    % Signal correlation
    plot([0 1],[Rsig Rsig],'-','Color',colors(4,:),'LineWidth',4)
    
    % Also correlate the MPH shapes for each stim
    if ist>1 && ist<7
        r = corrcoef(ThisData(i_sorted(Un1),1:ceil(1000/AMrates(ist-1)),ist-1),ThisData(i_sorted(Un2),1:ceil(1000/AMrates(ist-1)),ist-1));
        plot([0 1],[r(1,2) r(1,2)],'-','Color','b','LineWidth',4)
    elseif ist==1
        r = corrcoef(FR_Warn(i_sorted(Un1),:),FR_Warn(i_sorted(Un2),:));
        plot([0 1],[r(1,2) r(1,2)],'-','Color','b','LineWidth',4)
    end
    
    % Noise correlation
    plot(0.5, Rnoise_stim(ist) ,'x','Color',colors(4,:),'MarkerSize',20,'LineWidth',6)
    
    ylabel('Corr')
    set(gca,'xtick',[],'xlim',[0 1],'Color','none',...
        'ylim',[-1 1],'ytick',[-1 0 1])
    
    title(Info.stim_ID_key{ist})
end

suptitle(sprintf('Signal and noise correlations:  %i (%s) & %i (%s)',Clu1,UnI.SpkShape{Un1},Clu2,UnI.SpkShape{Un2} ))

savedir = fullfile(fn.figs,'NoiseCorr',SESSION);
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf,fullfile(savedir,sprintf('NC_clus_%i_%i',Clu1,Clu2)))


end