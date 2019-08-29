
fn = set_paths_directories;

% Load T
load(fullfile(fn.processed,'UnLUTcutDB.mat'),'-mat')

% Load Units files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Step through units in Table

for ic = 1:size(T,1)
    
    % Find corresponding iUn index
    iUn = find( strcmp(UnitInfo.Subject,T(ic,:).subject) & strcmp(UnitInfo.Session,T(ic,:).session) & [UnitInfo.Clu]==T(ic,:).clu );
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    
    % Load data files
    fprintf('Loading %s sess %s...\n',subject,session)
    clear TrialData Info Clusters Spikes
    filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    % Get spiketimes and shift based on calculated integration time 
    if exist('Spikes','var')                                 % >>> UMS <<<
        spiketimes = unique(round(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift));  %ms
    elseif exist('Clusters','var')                           % >>> KS  <<<
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(round(Clusters(iClu).spikeTimes * 1000 - spkshift)');
    end
    
    % Smooth FR 1 min window
    Stream_FRsmooth = convertSpiketimesToFR(spiketimes,length(SpoutStream),...
        TrialData.onset(1),TrialData.offset(1),60000,60000,'silence');
    
    
    % Find TD idx of first DB (or AC)
    Irrs = unique(TrialData.trID,'stable');
    Irrs = Irrs(Irrs>6);
    
    TDidx = find(TrialData.trID==Irrs(2),1,'first');
    
        
    % Plot smoothed FR to find time to cut off
    figure; plot(Stream_FRsmooth,'k')
    hold on
    plot(TrialData(TDidx,:).onset*[1 1],[0 max(Stream_FRsmooth)],'-g')
    
    keyboard
    continue
    
    %~~~~~~~~~~~~~~~
    % Optional code to run manually if FR clearly and abruptly decreases 
    
    Info.artifact(channel).trials = [Info.artifact(channel).trials; (TDidx:size(TrialData,1))'];
    save(fullfile(fn.processed,subject,sprintf('%s_sess-%s_Info',subject,session)),'Info','-v7.3');
    
    T.CutTrial(ic) = TDidx;
    save(fullfile(fn.processed,'UnLUTcutDB'),'T','-v7.3')
    %~~~~~~~~~~~~~~~
    
end


