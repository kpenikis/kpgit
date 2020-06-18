function plotPopClassEx
%
%  KP, 2020-06
%

% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
winlen       = 50;
convwin      = exp(-lambda*(1:winlen));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('shuffle')


%% Load data

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        rootdir = fullfile(fn.figs,'ClassAM');
        rawdata = 'CTTS_AM';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassSpeech');
        rawdata = 'CTTS_Speech_nonSim';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
end

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(rootdir,'RawData',rawdata)); %Cell_Time_Trial_Stim
CTTS = q.Cell_Time_Trial_Stim;


% CellTypes
iRS = find(UnitInfo.TroughPeak>0.43);
iNS = find(UnitInfo.TroughPeak<=0.43);



%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];

% Set figsavedir
figsavedir = fullfile(rootdir,whichStim,'PopExRasters');
if ~exist(figsavedir,'dir')
    mkdir(figsavedir)
end


%%

itr = 15;
ist = 2; 
ncells = 80;

t = 491:990;

hf=figure;
set(gcf,'Position',tallsmall)
hold on

CTTS(iRS(isnan(sum(CTTS(iRS,t,itr,ist),2))),:,itr,ist) = 0;
[y,x] = find(CTTS(iRS((1:ncells)+0),t,itr,ist));

% Raster
subplot(4,1,1:3)
plot([x x]',[-0.5*ones(size(y))+y 0.5*ones(size(y))+y]','k','LineWidth',2)
% plot(x',y','.k','MarkerSize',8)

xlim([0 500])
ylim([0 ncells])
set(gca,'tickdir','out','ticklength',[0.025 0.025],'ytick',[],'xtick',[],'Color','none')
box off

% PSTH
subplot(4,1,4)
plot(1:500,conv(sum(CTTS(iRS,t,itr,ist),1),convwin,'same'),'k','LineWidth',2)

xlim([0 500])
ylim([0 5])
set(gca,'tickdir','out','ticklength',[0.025 0.025],'Color','none')
box off

% plot(1:500,nansum(nansum(CTTS(iRS,501:1000,:,ist),3),1),'m','LineWidth',2)

title(sprintf('stim: %i, tr: %i',ist,itr))


print_eps_kp(hf,fullfile(figsavedir,sprintf('%s_Stim%i_Tr%i_nRS%i',whichStim,ist,itr,ncells)))

end




