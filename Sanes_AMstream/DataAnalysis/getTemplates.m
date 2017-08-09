function [T,ntUsed] = getTemplates(blocksN,minDur,ntUsed,spiketimes,spl,lpn,amd,ArtifactFlag,fs,SoundData,GW)

global mspad

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
T = cell(2,length(blocksN));

for ib=1:length(blocksN)
    
    % Choose random trial from range, but check that it hasn't been used
%     nt=0;
%     while ismember(nt,ntUsed(:,ib))
        nt = randi(blocksN(ib),1);
%     end
    % Once have a template trial, add it to matrix to track which used
    ntUsed(it,ib) = nt;
    
    
    % Get trial times
    [~,~,bkStart_ms,bkStop_ms] = get_blockOnsets( SoundData,ib,...
                                    spl,lpn,amd,ArtifactFlag,fs);
    
    % Get spiketimes during the selected trial, adjusted to start time 0
    try
    sp=[]; 
    sp = spiketimes( spiketimes>=(bkStart_ms(nt)-mspad) & spiketimes<(bkStop_ms(nt)+mspad) ) - bkStart_ms(nt) + 1 + mspad;
    
    % Convert from spiketimes to binary vector
    sp01 = zeros(1,bkStop_ms(nt)-bkStart_ms(nt)+2*mspad);
    sp01(sp) = 1;
    
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
    T{1,ib} = sp_conv;
    T{2,ib} = nt;
    
end

end