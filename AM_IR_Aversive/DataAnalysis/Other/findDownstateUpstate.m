function TRANS = findDownstateUpstate(SUBJECT, SESSION)
%
%
%  KP, 2018-03
%


%% Load data files

fn = set_paths_directories(SUBJECT,SESSION,1);

filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

[UnitInfo, UnitData, Info, TrialData, Clusters, ~, artifactTrs ] = collectRasterDataSession(SUBJECT,SESSION);


%% 

theseUnits = find((cellfun(@numel,{Clusters.spikeTimes})/length(SpoutStream)*1000) > 1);

PopRaster = zeros(numel(theseUnits),length(SpoutStream));
for iUn = theseUnits
    PopRaster(iUn==theseUnits,round(Clusters([Clusters.clusterID]==UnitInfo.Clu(iUn)).spikeTimes*1000)') = 1;
end
sumPop = sum(PopRaster,1);


%% 

tstep  = 1;
twin_D = 50;
twin_U = 100;
nspk_D = 2;
nspk_U = 12;

DOWN = zeros(1,length(sumPop));
UP = zeros(1,length(sumPop));
TRANS = zeros(1,length(sumPop));
t0 = twin_D+1;
while (t0+1+twin_U)<=length(sumPop)
    
    if (sum(sumPop((t0-twin_D):t0)) < nspk_D) 
        DOWN(t0) = 1;
    end
    if (sum(sumPop((t0+1):(1+t0+twin_U))) > nspk_U) 
        UP(t0+1) = 1;
    end
    if (sum(sumPop((t0-twin_D):t0)) < nspk_D) && (sum(sumPop((t0+1):(1+t0+twin_U))) > nspk_U)
        TRANS(t0) = 1;
    end
    t0 = t0+tstep;
    
end
sum(diff(TRANS)==1)


% figure;
% plot(DOWN,'b')
% hold on
% plot(UP,'r')
% plot(TRANS,'g')
% ylim([-0.2 1.2])

iTrans = find(diff(TRANS)==1);

figure; 
hold off
for it = 1:length(iTrans)
    imagesc(PopRaster(:,iTrans(it)+[-100:100]))
    colormap(bone)
    pause(1.5)
end


end %function




