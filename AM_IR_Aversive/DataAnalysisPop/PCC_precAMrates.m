
global AMrates

fn = set_paths_directories;

% Matched MP data
q = load(fullfile(fn.processed,'MatchedMPs','MatchedMPdata'));
Data = q.Data;
clear q

% Population d' results
load(fullfile(fn.processed,'MPcontext','Pop','Pop_dPrimes'))

% SU classifier results
% q = load(fullfile(savedir,'MPcontextSUclassRes'));
% Res = q.Res;
% clear q


% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [  ...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;

% Plot discriminability as a function of distance of rates preceding Irr MP
% from itself
 

% Quantify preceding AM rates

iUn = 0;
np  = 0;
while np<4
    iUn=iUn+1;
    np = size(Data(iUn,2).data,1);
end

Prev500 = nan(4,5);
Prev100 = nan(4,5);

for irate = 1:5
    for ip = 1:np
        Prev500(ip,irate) = log2(Data(iUn,irate).data(ip,2).PrevAMrt500) - log2(AMrates(irate));
        Prev100(ip,irate) = log2(Data(iUn,irate).data(ip,2).PrevAMrt100) - log2(AMrates(irate));
    end
end


hf=figure;
set(gcf,'Position',fullscreen)

subplot(2,2,1)
plot(abs(Prev500),Pop_dPrimes,'ok','LineWidth',2)
axis square
ylabel('dprime')
xlabel('log distance of avg AM rate in preceding 500 ms')
xlim([0 5])

[r,p]=corr(abs(Prev500(:)),Pop_dPrimes(:));
title(sprintf('r=%0.2f, p=%0.2f',r,p))


subplot(2,2,2)
plot(abs(Prev100),Pop_dPrimes,'ok','LineWidth',2)
axis square
ylabel('dprime')
xlabel('log distance of avg AM rate in preceding 100 ms')
xlim([0 5])

[r,p]=corr(abs(Prev100(:)),Pop_dPrimes(:));
title(sprintf('r=%0.2f, p=%0.2f',r,p))


subplot(2,2,3)
hold on
for ir=1:5
    [prt,idx] = sort(abs(Prev500(:,ir)));
    plot(prt,Pop_dPrimes(idx,ir),'-','Color',colors(ir,:),'LineWidth',2)
end
axis square
ylabel('dprime')
xlabel('log distance of avg AM rate in preceding 500 ms')
xlim([0 5])

subplot(2,2,4)
hold on
for ir=1:5
    [prt,idx] = sort(abs(Prev100(:,ir)));
    plot(prt,Pop_dPrimes(idx,ir),'-','Color',colors(ir,:),'LineWidth',2)
end
axis square
ylabel('dprime')
xlabel('log distance of avg AM rate in preceding 100 ms')
xlim([0 5])


savedir = fullfile(fn.processed,'MPcontext','Pop');
print_eps_kp(hf,fullfile(savedir,'MPcontext_dpPop_precAMrate'))

% 
% 
% figure; 
% subplot(2,2,1)
% plot((Prev500),Pop_dPrimes,'.k')
% axis square
% title('Prior 500 ms')
% ylabel('dprime')
% xlabel('mean AM rate')
% 
% subplot(2,2,2)
% plot((Prev100),Pop_dPrimes,'.k')
% axis square
% title('Prior 100 ms')
% ylabel('dprime')
% xlabel('mean AM rate')
% 
% subplot(2,2,3)
% plot((Prev500),Pop_dPrimes,'-')
% axis square
% title('Prior 500 ms')
% ylabel('dprime')
% xlabel('mean AM rate')
% 
% subplot(2,2,4)
% plot((Prev100),Pop_dPrimes,'-')
% axis square
% title('Prior 100 ms')
% ylabel('dprime')
% xlabel('mean AM rate')
% 



