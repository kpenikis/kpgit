

% X: Make population matrix of binned spikecounts --  time(bin) x units 

% for Y = each bin vector 
% calculate malanobis distance of each bin observation from the rest
% d2 = mahal(Y,X)

% What does the distribution of d2 values?
% Pull the larger d2 values and ask what the stimulus was at that time


fn = set_paths_directories('','',1);

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitInfo = q.UnitInfo;
UnitData = q.UnitData;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

SESSIONS = unique(UnitInfo(:,1:2));

for ii = 2
    
    subject   = SESSIONS{ii,1}{:};
    session   = SESSIONS{ii,2}{:};
    
    % Load Info and TrialData files
    clear Info TrialData
    filename = sprintf( '%s_sess-%s_Info'      ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); load(fullfile(fn.processed,subject,filename));
    clear Clusters Spikes
    filename = sprintf( '%s_sess-%s_Spikes'    ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    if ~exist('Clusters','var') || numel(Clusters)<5
        continue
    end
    
    FullSpikeMatrix = zeros(numel(Clusters),round(max(vertcat(Clusters.spikeTimes)*1000)));
    
    for iClu = 1:numel(Clusters)
        
        % Get spiketimes (KS)
        spiketimes = unique(round(Clusters(iClu).spikeTimes*1000 - spkshift)');
        
        % Make gaussian smoothed version of activity
        Stream_FRsmooth = convertSpiketimesToFR(spiketimes,...
            length(FullSpikeMatrix),TrialData.onset(1),TrialData.offset(1),100,100,'silence');
        
        % Add data to matrix
%         FullSpikeMatrix(iClu,spiketimes) = 1;
        FullSpikeMatrix(iClu,:)          = Stream_FRsmooth;
        
    end
    
    ms1 = TrialData.onset(1);
    FullSpikeMatrix = FullSpikeMatrix(:,ms1:end);
    
    
    mdists = mahal(FullSpikeMatrix',FullSpikeMatrix');
    
    figure; 
    histogram(mdists,'DisplayStyle','stairs','Normalization','cdf')
    xlim([0 100])
    
    
    figure; 
    set(gcf,'DefaultAxesColorOrder',cmocean('algae'))
    set(gca,'NextPlot','replacechildren')
    for iex = 0.95:-0.05:0.5
        
        [pks,locs] = findpeaks(mdists,'MinPeakHeight',quantile(mdists,iex));
        locs = locs(1:end-1);
        
        foo = SoundStream(locs+ms1-1+[-100:2:200]);
        plot(median(foo,1)','LineWidth',2)
        
    end
    
    figure; 
    plot(mdists,'k')
    hold on
    plot(SoundStream(ms1:end).*quantile(mdists,0.95))
    
end











