function Result = getLatDur(theseCells,ist,CTTS,Result,CReach)
% Result = getLatDur(theseCells,ist,CTTS,Result,CReach)
% Called by LatDur_CTTS
% KP, 2020-03

global dpSU AnWin whichLat FRdiffThresh


% Edit Results struct
% all latency options
% try half width of findpeaks fct
% stack Spec with Rest in histogram 


for iu = 1:numel(theseCells)
    
    pk_time   = nan;
    t_onset   = nan;
    t_offset  = nan;
    t_mid     = nan;
    t_up      = nan;
    
    % Get latency and duration
    Data = permute(mean(CTTS(theseCells(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
    
    % Find half max FR and peak events above it
    peakThresh = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /4) + min(Data(:,AnWin));
    if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
        continue
    end
    [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',peakThresh);
    
    if isempty(PKS)
        continue
    end
    
    % Remove peaks too early
    %             PKS(LOCS<5) =[];
    %             LOCS(LOCS<5)=[];
    
    % ASSUME just one peak
    [~,ipk] = max(PKS);
    pk_height = PKS(ipk);
    pk_time   = LOCS(ipk);
    
    
    % Find time just before surpassing halfMax
    %         for ims = [pk_time:-1:1 pdms:-1:pk_time]
    for ims = ( AnWin(1) - 1 +  pk_time ) : -1 : 1
        if Data(:,ims)<peakThresh
            t_onset = ims - AnWin(1) - 1 ;
            break
        end
    end
    
    % Find time just after falling below halfMax
    %         for ims = [pk_time:pdms 1:pk_time]
    for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2)
        if Data(:,ims)<peakThresh
            t_offset = ims - AnWin(1) - 1 ;
            break
        end
    end
    
    % Calculate peak duration
    if t_offset>=t_onset
        t_up = t_offset-t_onset;
    end
    
    % Time of midpoint of peak
    t_mid = t_onset+t_up/2;
    
    % Save data for statistics
    Result(theseCells(iu),1,ist)  = 1;
    
    switch whichLat
        case 'pk_time'
            thisLat = pk_time; %t_onset; %t_mid;
        case 'onset'
            thisLat = t_onset; %t_mid;
    end
    
    Result(theseCells(iu),2,ist)  = thisLat;
    Result(theseCells(iu),3,ist)  = t_up;
    Result(theseCells(iu),4,ist)  = ceil(pk_height);
    Result(theseCells(iu),5,ist)  = CReach.dprime(theseCells(iu)); %mean/task d'
    Result(theseCells(iu),6,ist)  = dpSU(theseCells(iu),ist); %stimulus d'
    
end %iu

end