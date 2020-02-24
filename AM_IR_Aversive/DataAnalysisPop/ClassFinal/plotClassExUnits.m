
close all
%%%%%%%%%%%%%%%%%%%%%
whichClass  = 'Full';
whichStim   = 'Speech';
%%%%%%%%%%%%%%%%%%%%%
crAmt = 0.001;
%%%%%%%%%%%%%%%%%%%%%
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%

% Data settings
fn = set_paths_directories('','',1);

savedir = fullfile(fn.figs,'ClassResults',whichClass);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

%%%%%%%%%%%%%%%%%%%%%
switch whichStim
    case 'AC'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassAM','AC',whichClass);
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
        % Load SU data
        q=load(fullfile(datadir,'each','CR_each.mat'));
        CReach = q.CR;
        clear q
        
        % Also get CTTS
        [CTTS,theseCells] = recallDataParams('AC','each',12,1);
        
        
    case 'Speech'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassSpeech','Speech',whichClass);
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
        % Load SU data
        q=load(fullfile(datadir,'each','CR_each.mat'));
        CReach = q.CR;
        clear q
        
        % Also get CTTS
        [CTTS,theseCells] = recallDataParams('Speech','each',12,1);
        
end


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%% ~~~~ d'/PC vals for each stimulus

% CellTypes
iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);


[dps,iSUdps] = sort(CReach(iRS,:).dprime,'descend');

% UnitData(theseCells(iRS(iSUdps(6))))

thisUn = 23;
dps(thisUn)
unitstring = sprintf('%s_%i',UnitData(theseCells(iRS(iSUdps(thisUn)))).Session,UnitData(theseCells(iRS(iSUdps(thisUn)))).Clu);

%%%%%%%%%%%%%%%%%%%%%
% 	RASTERS
%%%%%%%%%%%%%%%%%%%%%
hfr=figure; 
suptitle(unitstring)
set(hfr,'Position',widescreen)

for ist = 1:size(CTTS,4)
    
    subplot(1,size(CTTS,4),ist)
    hold on
    
    for it=1:20
        sp = find( CTTS(iRS(iSUdps(thisUn)),:,it,ist) );
        plot([sp; sp],it*ones(size([sp; sp]))+[-0.5; 0.5],'-','Color',[1 1 1]*0.2,'LineWidth',2)
        set(gca,'Color','none','ytick',[])
    end
    ylim([0 it+1.5])
end


print_eps_kp(hfr,fullfile(savedir,['ExUnit_' unitstring '_rasters']))


%%%%%%%%%%%%%%%%%%%%%
%    PSTH
%%%%%%%%%%%%%%%%%%%%%
[CTTS,theseCells] = recallDataParams('Speech','each',12,convwin);

hfp=figure; 
suptitle(unitstring)
set(hfp,'Position',widescreen)

for ist = 1:size(CTTS,4)
    
    subplot(4,size(CTTS,4),ist)
    
    plot_psth = mean(CTTS(iRS(iSUdps(thisUn)),:,:,ist),3,'omitnan')*1000;
    
    plot(plot_psth,'-','Color',[1 1 1]*0.2,'LineWidth',2)
    set(gca,'Color','none','ytick',[],'xtick',[])
    box off
    ylim([0 100])
end

print_eps_kp(hfp,fullfile(savedir,['ExUnit_' unitstring '_psths']))


% Conf Mat

ConfMat = mean(CReach(iRS(iSUdps(thisUn)),:).Results{:},3);
muPC    = mean(diag(ConfMat))*100;
dprime  = norminv(mean(diag(ConfMat)),0,1) - norminv(mean(ConfMat(~diag(diag(ConfMat)))),0,1);

nStim = size(ConfMat,1);

% Switch order of stimuli
theseStim  = [4 3 2 1 5 6 7 8]; 
ConfMat = ConfMat(theseStim,theseStim);

ConfMat(1:nStim+1:nStim*nStim) = -1.*ConfMat(1:nStim+1:nStim*nStim);

% Plot
hfc = figure;
imagesc(ConfMat)
axis square
caxis([-1 1])
cmocean('curl','pivot',0)
colorbar
ylabel('True stim')
xlabel('Assigned')
set(gca,'xtick',[],'ytick',[])
box off
axis off
title(sprintf('%s: %0.1f%%, d''=%0.2f',unitstring,muPC,dprime))

print_eps_kp(hfc,fullfile(savedir,['ExUnit_' unitstring '_confmat']))




keyboard


