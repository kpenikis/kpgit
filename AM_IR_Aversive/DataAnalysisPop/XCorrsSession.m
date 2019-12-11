function XCorrsSession(SUBJECT,SESSION)

fn = set_paths_directories('','',1);
filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

% Load unit files
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
% Filter to just this session
iUnOrig  = find(strcmp(UnitInfo.Session,SESSION) & strcmp(UnitInfo.Subject,SUBJECT));
UnitData = UnitData(iUnOrig);
UnitInfo = UnitInfo(iUnOrig,:);

close all

for i1 = 1:numel(Clusters)-1
    for i2 = (i1+1):numel(Clusters)
        
        [c,lags] = xcorr(round(Clusters(i1).spikeTimes*1000)',round(Clusters(i2).spikeTimes*1000)',100);
        
        pk=[];
        pk=findpeaks(c);
        
        if ~isempty(pk)
            figure;
            plot(lags,c)
            title(sprintf('Clus %i & %i',i1,i2))
            
            
            u1 = find(UnitInfo.Clu==Clusters(i1).clusterID);
            u2 = find(UnitInfo.Clu==Clusters(i2).clusterID);
            
            UnitInfo([u1 u2],:)
            
            % Plot PSTHs for an Irr stim
            
            
            
        end
        
    end
end



end