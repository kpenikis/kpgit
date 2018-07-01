function [score,rel] = calculateCorrScore(dat,winLen,stimstr,drawPLOTS)

fs=1e3;
% dat should be a nxm matrix representing the binary spike times. 
% n is the number of trials 
% m is time bins, you can bin this if youd like and put in 4's etc. code
% will still work, but the tick marks on figure 1 will no longer represent
% spikes accurately
% shows spikes. 
n=size(dat,1);
timeVec=(1:size(dat,2))/fs*1e3;

%% Define your Smoothing Bins
% winLen=30;%this is relative to bins. i.e. if your samping rate is 1 kHz, 30 is a 30 ms gaussian
%this is the most important value!!! set to something that makes sense.

%gaussian window
window=gausswin(winLen);
window=window-min(window);
window=window/sum(window);

%rectangular window
% window=rectwin(winLen);
% window=window/sum(window);



%% Calculate score

miss=0;
s=nan(size(dat));%intialize smoothed matrix
thresh=Inf;%we were being biased by large values. put a ceiling them so that they dont drive the corr score. Inf means no threshold
for ii=1:n
    sig=conv(dat(ii,:),window,'same')*1e3;%MAKE SMOOTHIE THINGIE
    sig(sig>thresh)=thresh;%put a ceiling on large bursting events
    s(ii,:)=sig;%store
    if ~any(dat(ii,:))
        miss = miss+1;
    end
end

% Do Correlations
c=nan(n);
for ii=1:n% for all trials
    for j=1:ii-1%up until that trial
        c(ii,j)=corr(s(ii,:)',s(j,:)');%CORRELATE SMOOTHIE THINGIE
    end
end

score=nanmedian(c(~isnan(c)));

rel = (n-miss)/n;



if drawPLOTS
    
%% Plot Data and smoothed firing rate.

figure;
subplot(2,2,[1,3]);hold on;
o=1;%offset for plotting
scale=40;%240-2*winLen;%scale for plotting. play with these numbers

for ii=1:n
    plot(timeVec,n-(s(ii,:)/scale+o*ii-1),'color','k')%plot smoothed correlation
    %the "n-" trick is just to get the first trial on top. 
    %put tick marks
    spikes=find(dat(ii,:));
    for spk=1:length(spikes)
        line(timeVec(spikes(spk))*[1,1],n-([0,.5]+o*ii-1),'color',.7*[1,1,1])
    end
end
title('Spike Raster & Firing Rate')
xlabel('time (ms)')
ylabel('trial #')
ylim([0 o*n])
xlim([0,size(dat,2)]);
set(gca,'ydir','reverse')


% Show correlations
subplot(2,2,2)
imagesc(c)
colorbar;
% Set colormap so NaNs show up as black.
cmap=bone(256);
cmap(1,:)=[0,0,0];
colormap(cmap);
set(gca,'clim',[-1,1])%normalize color axis. you can take this out
title('X Corrs')
ylabel('trial #')
xlabel('trial #')
axis square
set(gca,'ydir','normal')


% Plot correlations and get score
subplot(2,2,4)%get score
histogram(c(~isnan(c)),-1:.1:1)%show all corellations
axis tight;
xlabel('r')
ylabel('N Pairs')

line(score*[1,1],ylim)
title(['Score = ' num2str(score,2)])

suptitle(stimstr)


end


end