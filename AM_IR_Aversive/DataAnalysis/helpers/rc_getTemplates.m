function [T,ntUsed] = rc_getTemplates(StreamSpikes,GW,AllStimStarttimes,minDur,ntUsed,mspad)


% Choose one template trial for each block
%  T{1,ib} = conv spiketimes
%  T{2,ib} = nt

% Add empty row to ntUsed
it = size(ntUsed,1)+1;
ntUsed(it,:) = zeros(1,size(ntUsed,2));

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
        nt = randi(length(Starttimes),1);
%     end
    % Once have a template trial, add it to matrix to track which used
    ntUsed(it,is) = nt;
    
    
    % Get trial times
    bkStart_ms = Starttimes(nt);
    bkStop_ms = bkStart_ms+minDur-1;
    
    % Get spiketimes during the selected trial, adjusted to start time 0
    sp01=[]; 
    sp01 = StreamSpikes((bkStart_ms-mspad):(bkStop_ms+mspad));
    try
%     sp=[]; 
%     sp = spiketimes( spiketimes>=(bkStart_ms(nt)-mspad) & spiketimes<(bkStop_ms(nt)+mspad) ) - bkStart_ms(nt) + 1 + mspad;
%     
%     % Convert from spiketimes to binary vector
%     sp01 = zeros(1,bkStop_ms(nt)-bkStart_ms(nt)+2*mspad);
%     sp01(sp) = 1;
    
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