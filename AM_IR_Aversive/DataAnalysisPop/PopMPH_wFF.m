function PopMPH_wFF
% 
% PopMPH
%
% 
%  Intended to help categorize response types.
% 

global AMrates useFR Boundaries


%~~~~~~~~~~~~~~~~~~~~~
useFR    =   'log'; 
Rotate   =    1;
%~~~~~~~~~~~~~~~~~~~~~
clipZ    =    0;
%~~~~~~~~~~~~~~~~~~~~~

close all

%% Load data

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load PopMPH data
savedir = fullfile(fn.figs,'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
load(fullfile(savedir,'MPHdata.mat'))

if size(FR_vec,1) ~= size(UnitInfo,1)
    warning('saved data files have different n units')
    plotPopNormMPH(1)
    load(fullfile(savedir,'MPHdata.mat'))
end

% Get FF data
FanoFs = FanoFactorDistributions(0);


%% Data settings

switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5 1.75]) ];
end


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',8)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [ 250 250 250;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];


%% Filter data to just one session

% keyboard
% SUBJ = 'AAB_265054';
% SESS = 'Apr02-AM';
% 
% theseUnits = strcmp(UnitInfo.Subject,SUBJ) & strcmp(UnitInfo.Session,SESS);
% 
% ThisData = ThisData(theseUnits,:,:);
% FR_Warn  = FR_Warn(theseUnits,:);
% 
% FanoFs   = FanoFs(theseUnits,:);


%%

% How do you want to sort the units?

[sortdata,i_sorted] = sort(FanoFs(:,1));


% Warn stim

hfw = figure;
set(hfw,'Position',fullscreen)
hold on

subplot(1,5,1);
plotWarnPopH


% Also plot FF for this stim and for UNMOD
subplot(1,5,2:3);
plot(FanoFs(i_sorted(1:ndp),1),1:ndp,'k')
set(gca,'ydir','reverse','Color','none')
ylim([0.5 ndp+0.5])
% set(gca,'ytick',r,'yticklabel',Boundaries)
box off
axis fill
title('FF during UNMOD')
xlabel('FF')
xlim([0 4])


% Save figure
savename = sprintf('PopMPH_%s_Unmod_sortFF',useFR);
if clipZ>0
    savename = [savename '-clipZ' num2str(10*clipZ)];
end

print_eps_kp(hfw,fullfile(fn.figs,'FF',savename))
set(hfw,'PaperOrientation','landscape')
print(hfw,'-dpdf','-r500','-fillpage', fullfile(fn.figs,'FF',savename))




for ir=1:5
    
    hf(ir) = figure;
    set(gcf,'Position',fullscreen)
    hold on

    subplot(1,5,1);
    plotThisMPH
    
    
    % Also plot FF for this stim and for UNMOD
    subplot(1,5,2:3);
    plot(FanoFs(i_sorted(1:ndp),1),1:ndp,'k')
    set(gca,'ydir','reverse','Color','none')
    ylim([0.5 ndp+0.5])
%     set(gca,'ytick',r,'yticklabel',Boundaries)
    box off
    axis fill
    title('FF during UNMOD')
    xlabel('FF')
    xlim([0 4])
    
    subplot(1,5,4:5);
    plot(FanoFs(i_sorted(1:ndp),ir+1),1:ndp,'k')
    set(gca,'ydir','reverse','Color','none')
    ylim([0.5 ndp+0.5])
%     set(gca,'ytick',r,'yticklabel',Boundaries)
    box off
    axis fill
    title(['FF during ' num2str(AMrates(ir)) ' Hz'])
    xlabel('FF')
    xlim([0 4])
    
    % Save figure
    savename = sprintf('PopMPH_%s_%ihz_sortFF',useFR,AMrates(ir));
    if clipZ>0
        savename = [savename '-clipZ' num2str(10*clipZ)];
    end
    
    print_eps_kp(gcf,fullfile(fn.figs,'FF',savename))
    set(gcf,'PaperOrientation','landscape')
    print(gcf,'-dpdf','-r500','-fillpage', fullfile(fn.figs,'FF',savename))
    
end %ir





end