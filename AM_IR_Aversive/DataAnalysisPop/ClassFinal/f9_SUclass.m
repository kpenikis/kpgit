function f9_SUclass
% Results of SU classification for manuscript. 
% Created from select plots from plots_0323.
%
% KP, 2020-05



% Clustered matrix of SU d' by stimulus
% Pairwise correlation of stimulus SU d' vectors
%   subsequently added:
% Average confusion matrices
% Boxplots of SU d' by stimulus

% close all

whichClass   = 'Full';
whichStim    = 'AC';


%% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',14)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];


%% Prepare data 

fn = set_paths_directories;

savedir = fullfile(fn.figs,'ClassResults',whichClass,'Fig9',whichStim);
if ~exist(savedir,'dir')
    mkdir(savedir)
end

switch whichStim
    case 'AC'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ AM ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassAM',whichStim,whichClass);
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        exSU   = 275;   %2A
%         exSU   = 133;   %2B
        
        theseStim  = 1:8; 
        istlab = {'Unmod' '2hz' '4hz' '8hz' '16hz' '32hz' 'IrrA' 'IrrB'};
        
    case 'Speech'
        %~~~~~~~~~~~~~~~~~~~~~~~~~ Speech ~~~~~~~~~~~~~~~~~~~~~~~~~
        
        datadir = fullfile(fn.figs,'ClassSpeech',whichStim,whichClass);
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        UnitData = q.UnitData;
        clear q
        
        exSU   = 32;
        
        theseStim  = [4 3 2 1 5 6 7 8]; 
        istlab = {'I''m' 'blab' 'I can''t' 'ber'  'which you' 'to be chop' 'please' 'trees'};
end

% Load SU data 
q=load(fullfile(datadir,'each','CR_each.mat'));
CReach = q.CR;
clear q


% Also get CTTS
[CTTS,theseCells] = recallDataParams(whichStim,'each',12);


iRS = find(UnitInfo(theseCells,:).TroughPeak>0.43);
% iNS = find(UnitInfo(theseCells,:).TroughPeak<=0.43);

% Example unit index
exSU = theseCells(iRS)==exSU;


%% d' by stimulus

dpmaxval = 3;
nStim    = size(CReach(1,:).Results{:},1);


% For each cell, get its PC & d' for each stimulus
pcStim = nan(length(iRS),nStim);
dpStim = nan(length(iRS),nStim);
Sparss = nan(length(iRS),1);
for ii = 1:length(iRS)
        
    ConfMat = mean(CReach(iRS(ii),:).Results{:},3,'omitnan');
    
    pcStim(ii,:) = diag(ConfMat)';
    dpStim(ii,:) = dp_from_ConfMat(ConfMat,0.01);
    Sparss(ii) = calculateSparseness(dpStim(ii,:)');
end


%% Population 

%==========================================================================
%                               SU d' 
%                    cluster cells by tuning prefs
%==========================================================================

DATA = round(dpStim(:,theseStim),1);
idplab = sprintfc('%d',1:size(DATA,1));

hf_foo=figure; 
set(hf_foo,'Position',widescreen)

% Calulate linkages
Y = pdist(DATA);     % cells X time
Z = linkage(Y,'ward');

leafOrder = fliplr(optimalleaforder(Z,Y));


% subplot(2,1,1)     % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'Orientation','top','Labels',idplab,'reorder',leafOrder);%,'Labels',idplab);
% set(gca,'tickdir','out')

% Plot matrix only clustered by tuning curve similarity
% subplot(2,1,2)
% imagesc( DATA(fliplr(outperm),:)' )
% caxis([-0.2 dpmaxval])

DATA_sc = DATA(fliplr(outperm),:)';
sortCells = fliplr(outperm);
close(hf_foo)

%%

DATA2 = DATA_sc;

% Calulate linkages
Y = pdist(DATA2);     % cells X time
Z = linkage(Y,'ward');

leafOrder = fliplr(optimalleaforder(Z,Y));


% Plot
hf_mat=figure; 
set(hf_mat,'Position',widescreen)

subplot(7,4,[1 5])     % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'Orientation','top','Labels',istlab,...
    'reorder',leafOrder);  %,'reorder',leafOrder);%,'Labels',idplab);
set(gca,'tickdir','out','Color','none')
xtickangle(45)


% Clustered by tuning curve similarity
subplot(7,4,(4*3+1):(4*6))
imagesc( DATA2(fliplr(outperm),:) )
cmocean('thermal')
caxis([0 dpmaxval])
set(gca,'ytick',1:numel(istlab),'yticklabel',istlab(fliplr(outperm)),...
    'Color','none','xtick',round(linspace(0,numel(iRS),6)))
xlabel('RS cells')
box off
hcb=colorbar;
box off
hcb.Label.String = 'd''';



% print_eps_kp(hf_mat,fullfile(savedir,'SUdp_stim_sorted'))


sortStim = fliplr(outperm);




%==========================================================================
%%                        SU d' vs n spikes
%           Cumulative fraction of spikes -- cell/stim pairs
%==========================================================================


% Get Nsp_per_tr

Nsp_per_tr = nan(length(iRS),nStim);

for ii = 1:length(iRS)
    for ist=1:nStim
        clear y
        
        CSdata = permute(CTTS(iRS(ii),:,:,ist),[3 2 1]);
        [~,y] = find(CSdata==1);
              
        Nsp_per_tr(ii,ist) = numel(y) / sum(~isnan(CSdata(:,500)));
    end 
end


% Plot

[dpALLsort,iALLsort] = sort(dpStim(:),'ascend');
NspkALL = Nsp_per_tr(:);

figure; hold on
yyaxis left
plot(1:numel(NspkALL),cumsum(NspkALL(iALLsort))./sum(NspkALL(iALLsort)),'b','LineWidth',3)
% plot(1:numel(NspkALL),NspkALL(iALLsort),'b','LineWidth',1.5)
box off
set(gca,'Color','none')
yyaxis right
plot(1:numel(NspkALL),dpALLsort,'.k','MarkerSize',15)
xlim([1 numel(NspkALL)])
box off
set(gca,'Color','none')


skw_dp = skewness(dpStim(:));
skw_ns = skewness(NspkALL);

[r,p] = corrcoef(dpALLsort,NspkALL(iALLsort));

title(sprintf('r^2=%0.1f%%   y_dp=%0.2f   y_nspk=%0.2f',100*(r(1,2)^2),skw_dp,skw_ns))

% print_eps_kp(gcf,fullfile(savedir,'SU_cmNspks_CSpairs'))



% Double cumulative (Lorenz curve):

% figure; hold on
% yyaxis left
% plot(1:numel(NspkALL),cumsum(NspkALL(iALLsort))./sum(NspkALL(iALLsort)),'b','LineWidth',3)
% box off
% set(gca,'Color','none')
% yyaxis right
% plot(1:numel(NspkALL),cumsum(dpALLsort)./sum(dpALLsort),'.k','MarkerSize',15)
% xlim([1 numel(NspkALL)])
% box off
% set(gca,'Color','none')

% ref for gini coefficient when some values negative (like d')
% https://www.jstor.org/stable/2662589?seq=5#metadata_info_tab_contents



%% Example unit

% Switch order of stimuli
ConfMat = mean(CReach(iRS(exSU),:).Results{:},3);
ConfMat = ConfMat(theseStim(sortStim),theseStim(sortStim));

figure;

subplot(2,2,1)

imagesc(ConfMat)
axis square
caxis([0 1])
cmocean('-gray')
colorbar
ylabel('True stim')
xlabel('Assigned')
set(gca,'xtick',1:nStim,'xticklabel',sortStim,'ytick',1:nStim,'yticklabel',sortStim,...
    'Color','none','tickdir','out')

subplot(2,2,3)
plot(1:8,dpStim(exSU,theseStim(sortStim)),'g.-','MarkerSize',20)
set(gca,'xtick',1:nStim,'xticklabel',istlab(sortStim),...
    'Color','none','tickdir','out')
xlim([0.5 nStim+0.5])
ylim([-0.2 3])
% box off

print_eps_kp(gcf,fullfile(savedir,['SU_exUn-' num2str(theseCells(iRS(exSU)))]))


end

