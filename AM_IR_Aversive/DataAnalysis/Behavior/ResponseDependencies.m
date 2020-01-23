

SUBJECT = 'AAB_265054';
SESSION = 'Apr07-AM';


% CYCLE THROUGH SESSIONS, MASTER PLOT FOR EACH ANIMAL


fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));



%% Compare Hit Rates between Warns that follow 4 vs 32 Hz 

it_Warn = find(TrialData.trID==1);

% figure;
SpoutScore_4  = [];
SpoutScore_32 = [];
HitScore = nan(numel(it_Warn),2);
for it = 1:numel(it_Warn)
    
    t0 = TrialData.onset(it_Warn(it));
    t1 = TrialData.offset(it_Warn(it));
    
    OnAtBeg  = sum(SpoutStream(t0+(0:499)))>300;
    OffAtEnd = sum(SpoutStream(t1+(-499:0)))<200;
    EasyScore = numel(find(diff(SpoutStream(t0:t1))))==1;
    
    if OnAtBeg && OffAtEnd 
        HitScore(it,1) = 1;
        HitScore(it,2) = TrialData.trID(it_Warn(it)-1);
    elseif OnAtBeg && ~OffAtEnd
        HitScore(it,1) = 0;
        HitScore(it,2) = TrialData.trID(it_Warn(it)-1);
%     elseif ~OnAtBeg && ~OffAtEnd
%         continue
%     else
%         hold on
%         plot(SpoutStream(t0:t1))
% %         keyboard
%         HitScore(it,1) = 1;
%         HitScore(it,2) = TrialData.trID(it_Warn(it)-1);
    end
    
    if TrialData.trID(it_Warn(it)-1)==3
        SpoutScore_4  = [SpoutScore_4;  SpoutStream(t0+(0:1499))];
    elseif TrialData.trID(it_Warn(it)-1)==6
        SpoutScore_32 = [SpoutScore_32; SpoutStream(t0+(0:1499))];
    end
end

figure;
plot(mean(SpoutScore_4,1))
hold on
plot(mean(SpoutScore_32,1))
ylim([0 1])


sum(HitScore(HitScore(:,2)==3,1)==1) / sum(HitScore(:,2)==3)
sum(HitScore(HitScore(:,2)==6,1)==1) / sum(HitScore(:,2)==6)



%% Compare FA rate for 2 Hz periods:  
%   first pdc (sep 4 and 32 hz) vs last period of Irr AC


% Periodic
it_2Pdc = find(TrialData.trID==2);

SpoutScore_2_4  =[];
SpoutScore_2_32 =[];
for it = 1:numel(it_2Pdc)
    
    t0 = TrialData.onset(it_2Pdc(it));
    
    if TrialData.trID(it_2Pdc(it)-1)==3
        SpoutScore_2_4  = [SpoutScore_2_4;  SpoutStream(t0+(0:999))];
    elseif TrialData.trID(it_2Pdc(it)-1)==6
        SpoutScore_2_32 = [SpoutScore_2_32; SpoutStream(t0+(0:999))];
    end
end


t22=round(sum(1000./[2 4 4 8 8 16 16 32 32]));

% Irregular
it_2Irr = find(TrialData.trID==7);

SpoutScore_2_AC =[];
for it = 1:numel(it_2Irr)
    t0 = TrialData.onset(it_2Irr(it));
    SpoutScore_2_AC  = [SpoutScore_2_AC;  SpoutStream(t0+t22-1+(0:499))];
end


figure;
plot(mean(SpoutScore_2_4,1))
hold on
plot(mean(SpoutScore_2_32,1))
plot(mean(SpoutScore_2_AC,1),'LineWidth',2)
ylim([0 1])






