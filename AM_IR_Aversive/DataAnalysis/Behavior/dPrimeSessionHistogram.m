function SESSIONS = dPrimeSessionHistogram

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

SESSIONS = unique(UnitInfo(:,1:2));
dPrimes  = nan(size(SESSIONS,1),1);

for ii = 1:size(SESSIONS,1)
    
    subject   = SESSIONS{ii,1}{:};
    session   = SESSIONS{ii,2}{:};
    
    if strcmp(subject,'WWWlf_253395') && strcmp(session,'MA')
        dPrimes(ii) = nan;
        continue
    end
    
    % Load Info and TrialData files
    clear Info TrialData
    filename = sprintf( '%s_sess-%s_Info'      ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    % Load behavior data
    try
        load(fullfile(fn.raw,subject,[Info.epData_fn '_behavior.mat']));
    catch
        load(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '_behavior.mat']));
    end
    
    
    % Calculate behavioral performance
    
    GO = find([Data.TrialType]==0);
    
    G_hit  = bitget([Data(GO).ResponseCode],Info(1).Bits.hit);
    G_miss = bitget([Data(GO).ResponseCode],Info(1).Bits.miss);
    
    HitRate = sum(G_hit)/numel(G_hit);
    if HitRate == 1
        HitRate = 1 - 1/(2*numel(G_hit));
    elseif HitRate == 0
        HitRate = 1/(2*numel(G_hit));
    end
    
    
    NOGO = find([Data.TrialType]==1);
    
    NG_cr = bitget([Data(NOGO).ResponseCode],Info(1).Bits.cr);
    NG_fa = bitget([Data(NOGO).ResponseCode],Info(1).Bits.fa);
    
    FARate = sum(NG_fa)/numel(NG_fa);
    if FARate == 0
        FARate = 1/(2*numel(NG_fa));
    elseif FARate == 1
        FARate = 1 - 1/(2*numel(NG_fa));
    end
    
    
    %~~~~~  d prime  ~~~~~%
    
    zHit = norminv(HitRate);
    zFA = norminv([FARate]);
    
    dPrimes(ii) = zHit-zFA;
    
    
    if dPrimes(ii) < 1.5
        keyboard
        % if find another session with low d', load Info, make all trials
        % artifact, save Info, add to top of loop to skip
    end
    
    
end


SESSIONS.dprime = dPrimes;


% Plot
figure;
histogram(dPrimes,0:0.25:4)
xlabel('dprime of session')
ylabel('Count')


