function datOut = convolveGauss(dat,winLen)

if nargin<2
    winLen=30;%this is relative to bins. i.e. if your samping rate is 1 kHz, 30 is a 30 ms gaussian
    %this is the most important value!!! set to something that makes sense.
end

% Make gaussian window
window=gausswin(winLen);
window=window-min(window);
window=window/sum(window);

datOut=nan(size(dat));%intialize smoothed matrix

for i=1:size(dat,1)
    sig=conv(dat(i,:),window,'same');%MAKE SMOOTHIE THINGIE
    datOut(i,:)=sig;%store
end

end