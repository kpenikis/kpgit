function SESSIONS = dPrimeSessionHistogram

% Load Unit data files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

SESSIONS = unique(UnitInfo(:,1:2));
dPrimes  = nan(size(SESSIONS,1),1);
HitRate  = nan(size(SESSIONS,1),1);
FARate  = nan(size(SESSIONS,1),1);

for ii = 1:size(SESSIONS,1)
    
    subject   = SESSIONS{ii,1}{:};
    session   = SESSIONS{ii,2}{:};
    
%     if strcmp(subject,'WWWlf_253395') && strcmp(session,'MA')
%         dPrimes(ii) = nan;
%         continue
%     end
    
    % Load Info and TrialData files
    clear TrialData Info
    filename = sprintf( '%s_sess-%s_Info'      ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    % Load behavior data
    clear q Data
    try 
        q=load(fullfile(fn.raw,subject,[Info.epData_fn '_behavior.mat']));
    catch
        if strcmp(session,'Apr15-AM')
            continue
        else
            q=load(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '_behavior.mat']));
%             keyboard
        end
    end
      
    if ~exist('q','var')
        keyboard
    end
    
    DataBeh = q.Data;
    InfoBeh = q.Info;
    
    
    % Calculate behavioral performance
    
    GO = find([DataBeh.TrialType]==0);
    
    G_hit  = bitget([DataBeh(GO).ResponseCode],InfoBeh(1).Bits.hit);
    G_miss = bitget([DataBeh(GO).ResponseCode],InfoBeh(1).Bits.miss);
    
    HitRate(ii) = sum(G_hit)/numel(G_hit);
    if HitRate(ii) == 1
        HitRate(ii) = 1 - 1/(2*numel(G_hit));
    elseif HitRate(ii) == 0
        HitRate(ii) = 1/(2*numel(G_hit));
    end
    
    
    NOGO = find([DataBeh.TrialType]==1);
    
    NG_cr = bitget([DataBeh(NOGO).ResponseCode],InfoBeh(1).Bits.cr);
    NG_fa = bitget([DataBeh(NOGO).ResponseCode],InfoBeh(1).Bits.fa);
    
    FARate(ii) = sum(NG_fa)/numel(NG_fa);
    if FARate(ii) == 0
        FARate(ii) = 1/(2*numel(NG_fa));
    elseif FARate(ii) == 1
        FARate(ii) = 1 - 1/(2*numel(NG_fa));
    end
    
    
    %~~~~~  d prime  ~~~~~%
    
    zHit = norminv(HitRate(ii));
    zFA = norminv([FARate(ii)]);
    
    dPrimes(ii) = zHit-zFA;
    HitRate(ii) = HitRate(ii)*100;
    FARate(ii)  = FARate(ii)*100;
    
    
    if dPrimes(ii) < 1.5
%         keyboard
        % if find another session with low d', load Info, make all trials
        % artifact, save Info, add to top of loop to skip
    end
    
    
end

round(mean(dPrimes,'omitnan'))
round([min(dPrimes) max(dPrimes)])

round(mean(HitRate,'omitnan'))
round([min(HitRate) max(HitRate)])

round(mean(FARate,'omitnan'))
round([min(FARate) max(FARate)])



SESSIONS.dprime = dPrimes;


% Plot
figure;
histogram(dPrimes,0:0.25:4)
xlabel('dprime of session')
ylabel('Count')


