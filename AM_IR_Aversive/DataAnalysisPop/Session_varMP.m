function Session_varMP(SUBJECT,SESSION)
%
%  Session_varMP(SUBJECT, SESSION )
%    Should compare different data from different MP groups.
%     -within Pdc context: across various Pdc MPs
%     -within Irr context: each Irr MP separately?
% 
%  KP, 2019-09
%



%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load matched MP data
q = load(fullfile(fn.processed,'MatchedMPs','MatchedMPdata'));
Data = q.Data;
clear q

% Find indices of units from this session
idxSess = find(strcmp(UnitInfo.Subject,SUBJECT) & strcmp(UnitInfo.Session,SESSION))';


irate = 2;
nt=20;
Iterations=5000;

FF_mean = nan(numel(idxSess),2);
FF_std  = nan(numel(idxSess),2);

for iUn = idxSess
    
    vertcat(Data(iUn,irate).data(:,1).indices);
    rasterPdc = vertcat(Data(iUn,irate).data( size(Data(iUn,irate).data,1) + (-1:0) ,1).raster);
    rasterALL = vertcat(rasterPdc,vertcat(Data(iUn,irate).data(:,2).raster));
    
    % Bootstrap 20 trials within and across contexts, and calculate FF
    
    FF_Pdc = nan(Iterations,1);
    FF_ALL = nan(Iterations,1);
    for ii = 1:Iterations
        
        trs_Pdc = randperm(size(rasterPdc,1),nt);
        trs_ALL = randperm(size(rasterALL,1),nt);
        
        FF_Pdc(ii) = var(mean(rasterPdc(trs_Pdc,:),2)) / mean(mean(rasterPdc(trs_Pdc,:),2));
        FF_ALL(ii) = var(mean(rasterALL(trs_ALL,:),2)) / mean(mean(rasterALL(trs_ALL,:),2));
        
    end
    
    figure; 
    subplot(2,1,1)
    plotSpread([FF_Pdc FF_ALL],'showMM',1)
    title(['Clu ' num2str(UnitInfo.Clu(iUn))])
    
    subplot(2,1,2)
    plot(mean(rasterALL,1),'k','LineWidth',4)
    xlim([0 250])
    
    
    % Mean and SEM FF for this unit
    FF_mean(iUn==idxSess,:) = mean([FF_Pdc FF_ALL],1,'omitnan');
    FF_std(iUn==idxSess,:)  = std([FF_Pdc FF_ALL],1,'omitnan');
    
end


figure; 
plot([0 15],[0 15],'k')
hold on
plot(FF_mean(:,1),FF_mean(:,2),'ok')
axis square
xlabel('FF Pdc')
ylabel('FF ALL')


end