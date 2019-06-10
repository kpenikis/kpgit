
close all

% Get RateStream (fs = 1 kHz; 1 ms)
idxms = round( 1 : Info.fs_sound/1000 : length(epData.streams.rVrt.data(1,:)) );
RateStream = double(epData.streams.rVrt.data(1,idxms));


% Get Phase0 timestamps (fs = 1 kHz; 1 ms precision)
[~,AMMIN] = findpeaks(-SoundStream./max(SoundStream),'MinPeakProminence',0.4); %ms of RMS minima
Phase0 = [AMMIN; RateStream(AMMIN+5)];

% Remove Warn segments
WarnOnsets = TrialData([TrialData.trID]==1,:).onset';
for io = 1:numel(WarnOnsets)
    [~,ido] = min(abs(Phase0(1,:)-WarnOnsets(io)));
    if diff(Phase0(1,ido+(0:1))) > 1400
        Phase0(:,ido) = [];
    end
end


% Find outlier periods
pctDiff = abs( (Phase0(2,1:end-1) - round(1000./diff(Phase0(1,:))) ) ./ Phase0(2,1:end-1) );
CUTOFF = 0.25;

% Check outlier periods
for irate=[2 4 8 16 32]
    hf(irate) = figure; hold on
    idx = find(Phase0(2,:)==irate);
    for imin = idx
        if imin==length(Phase0), continue, end
        plot(SoundStream(Phase0(1,imin):Phase0(1,imin+1)),'k','LineWidth',2)
        xlim([1 10+1000/irate])
        ylim([0 1.5])
    end
end
for irate=[2 4 8 16 32]
    figure(hf(irate)); hold on
    idx = find(Phase0(2,:)==irate);
    for imin = idx
        if imin==length(Phase0), continue, end
        
        if pctDiff(imin) > CUTOFF
            plot(SoundStream(Phase0(1,imin):Phase0(1,imin+1)),'r','LineWidth',2)
        end
    end
end

% Remove outlier periods
for irate=[2 4 8 16 32]
    idx = find(Phase0(2,:)==irate);
    fprintf('   %i of %i %i Hz pds removed\n', sum(pctDiff(idx(1:end-1))>0.1),numel(idx),irate)
end

Phase0(:,pctDiff>CUTOFF) = [];
