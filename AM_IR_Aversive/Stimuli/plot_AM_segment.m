

subject = 'WWWf_253400';
session = 'MA';
BLOCK   =  39;


fn = set_paths_directories(subject,' ',1);

% Load this block epData datafile
block_str = sprintf('Block-%i.mat',BLOCK);
datafile = fullfile(fn.raw,subject,block_str);
if ~exist(datafile,'file')
    warning('raw epData file not found.')
    keyboard
end
fprintf(' loading data file %s...',datafile)
clear epData;
load(datafile,'-mat'); %loads data struct: epData
if isempty(epData)
    error('data file did not load correctly!')
else
    disp('done.')
end



%%

fs_sound   = epData.streams.rVrt.fs;
SoundData  = epData.streams.rVrt.data(2,:);
RateData   = epData.streams.rVrt.data(1,:);
DepthData  = epData.streams.rVrt.data(3,:);
% clear epData


sampBEG = round(5*60*fs_sound);%ms
sampEND = round(5.26*60*fs_sound);%ms


% Set up plot
AMrates = [2 4 8 16 32];
colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0 ]./255;
colors = [colors; 0.7.*bone(2)];

set(0,'DefaultAxesFontSize',24)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

hf=figure; hold on
set(gcf,'Position',figsize1)

% Plot rates
hs1=subplot(6,1,1:3); hold on
set(gca,'yscale','log','xlim',[sampBEG sampEND])
box off
plot(sampBEG:sampEND,RateData(sampBEG:sampEND),...
    'k','LineWidth',2)
set(gca,'xtick',[],'ytick',[2 4 8 16 32])
ylabel('AM rate (Hz)')

% Plot waveform
hs2=subplot(6,1,4:6); hold on
set(gca,'yscale','linear','xlim',[sampBEG sampEND])
box off
set(gca,'xtick',[],'ytick',[])

try
    
windowsize=100;
for ii = sampBEG:windowsize:sampEND
    
    subplot(hs1); hold on
    if DepthData(ii+windowsize/2)==0
        plot(ii+[0:windowsize],ones(size([0:windowsize])),'-',...
        'Color','k','LineWidth',8)
    else
        rates = RateData(ii+[0:windowsize]);
        if numel(unique(rates))==1
            plot(ii+[0:windowsize],RateData(ii+[0:windowsize]),'-',...
                'Color',colors(AMrates==RateData(ii+windowsize/2),:),'LineWidth',8)
        else
            Nrates = numel(unique(rates));
            iRateStart = [1 find(abs(diff(rates))>0)+1 windowsize+1];
            for irs = 1:(numel(iRateStart)-1)
                idx = ii+[iRateStart(irs):iRateStart(irs+1)-1]-1;
                plot(idx,RateData(idx),'-',...
                    'Color',colors(AMrates==RateData(idx(1)),:),'LineWidth',8)
            end
        end
    end
    hold off
    
    subplot(hs2); hold on
    plot(ii+[0:windowsize],SoundData(ii+[0:windowsize]),'-',...
        'Color',colors(AMrates==RateData(ii+windowsize/2),:))
    hold off
end


catch
    keyboard
end

keyboard

savedir = fn.stim;
savename = 'RMS_Rate_Segment';

print_eps_kp(hf,fullfile(savedir,savename));
print_svg_kp(hf,fullfile(savedir,savename));







