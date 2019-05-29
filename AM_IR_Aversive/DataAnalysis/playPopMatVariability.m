
SUBJECT = 'AAB_265058';
SESSION = 'Jan17-AM';


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Filter Unit files to just this session and sort by baseline FR
UnitData = UnitData(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitInfo = UnitInfo(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT),:);
[~, baseFRrank] = sort([UnitData.BaseFR]);
UnitInfo = UnitInfo(baseFRrank,:);
UnitData = UnitData(baseFRrank);


% Next try ranking units by FR in bin, converting to p(rank) and plotting
% average or spread of that
% Approaching likelyhood of activity vector 

% Also try something with Mahalinobis distance?


%%%%

figure; hold off

for ii = 16:-1:12
    
    FullSpikeMatrix = zeros(2,length(SoundStream));
    N=0;
    for iUn = ii:numel(UnitData)
        
        N=N+1;
        
        channel = UnitData(iUn).Channel(1);
        clu     = UnitData(iUn).Clu(1);
        
        % Get spiketimes (KS)
        spiketimes = unique(round(Clusters(([Clusters.maxChannel]==channel & [Clusters.clusterID]==clu)).spikeTimes*1000)');
        
        % Make gaussian smoothed version as well
        [Stream_FRsmooth,Stream_zscore,Stream_spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),20,100,'silence');
        
        % Add data to matrix
        FullSpikeMatrix(N,:) = Stream_FRsmooth;
        
    end
    
    
    Q1 = quantile(FullSpikeMatrix,0.25);
    Q3 = quantile(FullSpikeMatrix,0.75);
    
    FullQuartDispMatrix = (Q3-Q1)./(Q3+Q1);
    
    plot(FullQuartDispMatrix)
    ylim([-0.5 1.5])
    title([num2str(numel(UnitData)-ii+1) ' units'])
    pause(2)
    
end

hold on
plot(SoundStream./max(SoundStream),'k')


%%%%



FullCVMatrix = std(FullSpikeMatrix,1,'omitnan')./mean(FullSpikeMatrix,1,'omitnan');
FullFFMatrix = var(FullSpikeMatrix,1,'omitnan')./mean(FullSpikeMatrix,1,'omitnan');


figure; plot(FullCVMatrix)
hold on
plot(FullFFMatrix)





