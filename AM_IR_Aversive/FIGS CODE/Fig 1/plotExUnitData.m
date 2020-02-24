function plotExUnitData(SUBJECT,SESS_LABEL,Clu1,Clu2)
% 

% Load necessary data 
fn = set_paths_directories(SUBJECT,SESS_LABEL);

%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',12)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallfull = [1 scrsz(4) scrsz(3)/4 scrsz(4)];

% Set figsavedir
figsavedir = fullfile(fn.figs,'WaveformShapes','ExUnitWFs');
if ~exist(fullfile(figsavedir),'dir')
    mkdir(fullfile(figsavedir))
end


%%
load(fullfile( fn.processed,SUBJECT, sprintf( '%s_sess-%s_Info', SUBJECT,SESS_LABEL)));
load(fullfile( fn.processed,SUBJECT, sprintf('%s_sess-%s_Spikes',SUBJECT,SESS_LABEL)));

try
    load(fullfile(fn.sessdata,'sorting','KS_Output.mat'));
    clu_id     = readNPY(fullfile(fn.sessdata,'sorting','spike_clusters.npy'));
    spiketimes = readNPY(fullfile(fn.sessdata,'sorting','spike_times.npy')); % these are in samples, not seconds
catch
    keyboard
end


% Check root directory
fbinary_parts = strsplit(rez.ops.fbinary,'/');
if strcmp(fbinary_parts{4},'GDFS')
    rez.ops.root = fullfile(fn.processed,SUBJECT,SESS_LABEL,'sorting');
    rez.ops.fbinary = fullfile(fn.processed,SUBJECT,SESS_LABEL,'sorting',fbinary_parts{end});
    chanMapFile_parts = strsplit(rez.ops.chanMapFile,'/');
    rez.ops.chanMapFile = fullfile(fn.sorting,chanMapFile_parts{end-1},chanMapFile_parts{end});
end
load(rez.ops.chanMapFile,'chanMap','connected', 'xcoords', 'ycoords','kcoords');

fprintf('Fetching raw waveforms... ')
tic
gwfparams.dataDir = rez.ops.root;        % KiloSort/Phy output folder
% gwfparams.fileName = fbinary_parts{end}; % .dat file containing the raw    [(1:end-4) '_filt.dat'] ; 
fname = fbinary_parts{end}; % .dat file containing the raw    [(1:end-4) '_filt.dat'] ; 
gwfparams.fileName = [fname(1:strfind(fname,'.dat')-1) '_filt.dat'];
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = rez.ops.NchanTOT;        % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 500;                     % Number of waveforms per unit to pull out
gwfparams.spikeTimes = spiketimes;       % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = clu_id;        % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
gwfparams.theseUnits = [Clu1; Clu2];

wf = getWaveForms(gwfparams);
% wf.waveForms = [ Clu  x  iWaveform  x  Channel  x  sample ]
toc

idx1 = [Clusters.clusterID]==Clu1;
idx2 = [Clusters.clusterID]==Clu2;
Clusters(idx1).maxChannel
Clusters(idx2).maxChannel

ShankN = Clusters(idx1).shank;

switch ShankN
    case 1
        ChsToPlot = 0:11;
    case 2
        ChsToPlot = 12:23;
    case 3
        ChsToPlot = 24:35;
    case 4
        ChsToPlot = 36:47;
    case 5
        ChsToPlot = 48:59;
end

ChMap = 1+[0 11 1 10 2 9 3 8 4 7 5 6];

figure;
set(gcf,'Position',tallfull)
for ich = 1:numel(ChsToPlot)
    
    subplot(6,2,find(ich==ChMap))
    hold on
    
    theseWFs1 = permute(wf.waveForms(1,:,ich,:),[4 2 1 3]);
    theseWFs1 = theseWFs1(:,round(linspace(1,gwfparams.nWf,50)));
%     theseWFs1 = theseWFs1 - mean(theseWFs1(1:29,:));
    plot(theseWFs1,'Color',[0.5 0 0])
    
    theseWFs2 = permute(wf.waveForms(2,:,ich,:),[4 2 1 3]);
    theseWFs2 = theseWFs2(:,round(linspace(1,gwfparams.nWf,50)));
%     theseWFs2 = theseWFs2 - mean(theseWFs2(1:29,:));
    plot(theseWFs2,'Color',[0 0 0.5])
    
    plot(mean(permute(wf.waveForms(1,:,ich,:),[4 2 1 3]),2),'r','LineWidth',2)
    plot(mean(permute(wf.waveForms(2,:,ich,:),[4 2 1 3]),2),'Color',[0.4 0.5 1],'LineWidth',2)
    
    xlim([1 size(theseWFs2,1)])
    ylim([-12000 9000])
%     ylim([-10000 5000])
    set(gca,'Color','none','xtick',[])
    axis square
end

print_eps_kp(gcf,fullfile(figsavedir,sprintf('ExWFs_%s_%s_%i_%i_50', SUBJECT,SESS_LABEL,Clu1,Clu2)))

% WHAT ARE THE UNITS OF SPIKE AMPLITUDE?


%% 

% Autocorrelogram

xlimval = 200;

sp1   = round(Clusters(idx1).spikeTimes * 1000)';
sp2   = round(Clusters(idx2).spikeTimes * 1000)';

figure;
subplot(1,2,1);
histogram([-diff(sp1) diff(sp1)],[-xlimval:xlimval],'EdgeColor','none','FaceColor','r','FaceAlpha',1);
set(gca,'xlim',[-xlimval xlimval],'Color','none')
axis square
title(['Clu #' num2str(Clu1)])

subplot(1,2,2);
histogram([-diff(sp2) diff(sp2)],[-xlimval:xlimval],'EdgeColor','none','FaceColor',[0.4 0.5 1],'FaceAlpha',1);
set(gca,'xlim',[-xlimval xlimval],'Color','none')
axis square
title(['Clu #' num2str(Clu2)])

print_eps_kp(gcf,fullfile(figsavedir,sprintf('Autocorrs_%s_%s_%i_%i', SUBJECT,SESS_LABEL,Clu1,Clu2)))


% Amplitude over time
figure;
subplot(2,1,1);
hold on
plot(sp1,Clusters(idx1).Amplitudes,'.r');
plot(sp2,Clusters(idx2).Amplitudes,'.','Color',[0.4 0.5 1]);
ylim([0 40])

print_eps_kp(gcf,fullfile(figsavedir,sprintf('AmpTime_%s_%s_%i_%i', SUBJECT,SESS_LABEL,Clu1,Clu2)))


% PC space
iShank1 = find(chanMap(ChsToPlot+1)==Clusters(idx1).maxChannel);
iShank2 = find(chanMap(ChsToPlot+1)==Clusters(idx2).maxChannel);

iClu1 = clu_id==Clu1;
iClu2 = clu_id==Clu2;

% iShank1 = 10;
% iShank2 = 1;

figure;
plot(rez.cProjPC(iClu1,1,iShank1),rez.cProjPC(iClu1,1,iShank2),'.r')
hold on
plot(rez.cProjPC(iClu2,1,iShank1),rez.cProjPC(iClu2,1,iShank2),'.','Color',[0.4 0.5 1])
axis square
xlabel('Clu1 PC1')
ylabel('Clu2 PC1')
set(gca,'Color','none')

print_eps_kp(gcf,fullfile(figsavedir,sprintf('PCplot_%s_%s_%i_%i', SUBJECT,SESS_LABEL,Clu1,Clu2)))




end
