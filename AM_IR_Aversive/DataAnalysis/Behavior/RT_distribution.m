function RT_distribution


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units_250'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)
scrsz = get(0,'ScreenSize');    %[left bottom width height]

session='';

AllRTs  = nan(numel(UnitData),50);
AllRTs2 = nan(numel(UnitData),50);

for iUn = 1:numel(UnitData)
    
    if strcmp(session,UnitData(iUn).Session)
        continue
    end
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    
    % Load data files
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData 
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    
    % First for Warn
    idxWarn = find(TrialData.trID==1);
    
    RTsubj = [];
    for it = 1:numel(idxWarn)
        t0 = TrialData.onset(idxWarn(it));
        RTsubj = [RTsubj find(diff(SpoutStream(t0+[0:1500]))==-1,1,'first')];
        if find(diff(SpoutStream(t0+[0:1500]))==-1,1,'first')==0
            keyboard
        end
    end
    AllRTs(iUn,1:numel(RTsubj)) = RTsubj;
    
    % Also collect for 2 Hz
    idxStim = find(TrialData.trID==2);
    
    RTsubj = [];
    for it = 1:numel(idxStim)
        t0 = TrialData.onset(idxStim(it));
        RTsubj = [RTsubj find(diff(SpoutStream(t0+[0:1000]))==-1,1,'first')];
        if find(diff(SpoutStream(t0+[0:1000]))==-1,1,'first')==0
            keyboard
        end
    end
    AllRTs2(iUn,1:numel(RTsubj)) = RTsubj;
    
    
end %iUn


hf=figure;
histogram(AllRTs(~isnan(AllRTs)),0:10:1500,'FaceColor','k')
ylim([0 100])
text(20,95,[num2str(sum(sum(AllRTs<10))) '<10ms'])
xlabel('First time off spout')
ylabel('Count')
title('Warn')
print_eps_kp(hf,fullfile(fn.figs,'Behavior','RTdist_Warn'))

hf2=figure;
histogram(AllRTs2(~isnan(AllRTs2)),0:10:1000,'FaceColor','k')
ylim([0 10])
text(20,9.5,[num2str(sum(sum(AllRTs2<10))) '<10ms'])
xlabel('First time off spout')
ylabel('Count')
title('2 Hz')
print_eps_kp(hf2,fullfile(fn.figs,'Behavior','RTdist_2Hz'))



end