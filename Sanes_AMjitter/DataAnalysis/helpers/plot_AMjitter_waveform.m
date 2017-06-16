function plot_AMjitter_waveform(stimfn)

fn = set_paths_directories;
% stimfn = '0';
rV = jitter_LUT(4,stimfn);


fs = 1e3;
t  = cumsum([0 1./rV]);
s = round(t.*fs);
y  = [];
for ip = 1:numel(rV)
    tp=[]; tp = 0:(s(ip+1)-s(ip)-1);
    y = [y 1- 0.75*cos(2*pi*rV(ip)/fs*tp)];
end
%     hold on; plot(y)

% Remove first 1/4 of the first period
yclip = y(round(0.25/rV(1)*fs):end);

% Add sound during unmodulated portion
y_um = ones(1,round(212/1000*fs)); %0.8165.*
stim = [y_um yclip];

% Plot individual stim, for debugging
hF = figure;
scrsz = get(0,'ScreenSize');
set(hF,'Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2],...
    'Nextplot','add');
fill([1:length(stim) length(stim):-1:1] /fs*1000,[stim fliplr(-stim)],[0.3 0.3 0.3])
hold on
% xtx = 500:500:2000;
% for itx=xtx
%     plot([itx itx],[-2 2],'k--')
% end
ylim([-3 3])
xlim([0 length(stim)]/fs*1000)
xlabel('time (ms)')
title(stimfn)

set(hF,'PaperOrientation','landscape');
print(hF,fullfile(fn.stim,['waveform_' stimfn]),'-dpdf','-bestfit');


