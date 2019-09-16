function [T,ntUsed,TemplTrs] = rc_getTemplates(StreamSpikes,GW,AllStimStarttimes,minDur,ntUsed,mspad)

global exclOnset TempSize

% Choose one template trial for each block
%  T{1,ib} = conv spiketimes
%  T{2,ib} = nt

% Add empty row to ntUsed
it = size(ntUsed,1)+1;
ntUsed(it,:) = zeros(1,size(ntUsed,2));

TemplTrs = nan(TempSize,numel(AllStimStarttimes));

% Determine the half length (in samples) of the gaussian window. Use later
% for shifting convolved signal to align with correct time points.
lGausHalf = (length(GW)-1)/2;

% Make empty cell for output data
T = cell(2,numel(AllStimStarttimes));

for is=1:numel(AllStimStarttimes)
    
    Starttimes = AllStimStarttimes{1,is};
    if isempty(Starttimes)
        continue
    end
    % Choose random trial from range
    % only check that it hasn't been used if nTemplates is smaller than Ntrials
%     nt=0;
%     while ismember(nt,ntUsed(:,is))
        nt = randperm(length(Starttimes),TempSize);
%     end
    % Once have a template trial, add it to matrix to track which used
    ntUsed(it,is)  = nt(1);
    TemplTrs(:,is) = nt';
    
    % Get trial times
    bkStart_ms = Starttimes(nt);
    bkStop_ms = bkStart_ms+minDur-1;
    if exclOnset
        bkStart_ms = bkStart_ms+150;
    end
    
    % Get spiketimes during the selected trial, adjusted to start time 0
    sp01=nan(TempSize,length((bkStart_ms-mspad):(bkStop_ms+mspad))); 
    for ii = 1:TempSize
        sp01(ii,:) = StreamSpikes((bkStart_ms-mspad):(bkStop_ms+mspad));
    end
    sp01 = mean(sp01,1);
    
    try
    % Convolve with Gaussian window
    sp_conv = conv(sp01,GW,'full');
    
    % Correct for the time shift introduced by the convolution
    sp_conv = sp_conv(lGausHalf+1:end-lGausHalf);
    
    % Check by plotting convolved signal against spiketimes
%     figure; hold on
%     h = plot(sp_conv,'r-');
%     y = zeros(length(sp),1);
%     h2 = plot(sp,y,'r.');
    
    % Remove extra time added for padding AND cut block to min block dur
    sp_conv = sp_conv(mspad+(1:minDur));
    
    % Check by plotting again
%     figure; hold on
%     h = plot(sp_conv,'b-');
%     y = zeros(length(sp),1);
%     h2 = plot(sp-mspad,y,'b.');
%     xlim([1 minDur])
    
    catch
        keyboard
    end
    
    % Finally, add the data to the output array
    T{1,is} = sp_conv;
    T{2,is} = nt;
    
end

end