
% close all
useFR    =   'log'; 


% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


% Load RepSpeech templates
k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
kfns = fieldnames(k);

Duration = max(structfun(@length,k));
Segs     = 1:numel(structfun(@length,k));
SegDurs  = structfun(@length,k);


savedir = fullfile(fn.figs,'PopulationSpeechSegments');
load(fullfile(savedir,'RepSpeech'))

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end

%% Set up

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',8)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

%%

% Data settings
useFR    =   'z'; 

switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
%         Boundaries  = [-1 round(10.^[0.5]) ];
end

clipVal = 100;

iUn_GoodData = find(sum(sum(~isnan(ThisData),2),3)==sum(SegDurs));

GoodData = [];
for ist = 1:size(ThisData,3)
    
%     T = clusterdata(zFR_vec(:,:,ir),10);

    pdms = SegDurs(ist);
        
    % Remove empty units
%     iUn_GoodData = find(sum(~isnan(ThisData(:,:,ist)),2)==pdms);
    GoodData = [GoodData ThisData(iUn_GoodData,1:pdms,ist)];
    
end
GoodData(GoodData> clipVal)  =  clipVal;
GoodData(GoodData<-clipVal)  = -clipVal;


% Calulate linkages
Y = pdist(GoodData);     % cells X time
Z = linkage(Y,'ward');

leafOrder = fliplr(optimalleaforder(Z,Y));

% Plot dendrogram
figure;
set(gcf,'Position',fullscreen)

subplot(1,6,4:6);      % ,'ColorThreshold',40 ,'reorder',leafOrder
[hd,tvals,outperm] = dendrogram(Z,0,'Orientation','right','Labels',cellfun(@num2str, num2cell(iUn_GoodData),'UniformOutput',false));
set(gca,'tickdir','out')

RespCluLabels = cluster(Z,'maxclust',10);

% Plot MPH
subplot(1,6,1:3);

imagesc( GoodData(fliplr(outperm),:) )

caxis([-0.5 4])
cmocean('balance','pivot',0)
xlim([0 size(GoodData,2)])
ylim([0 size(GoodData,1)+1])
set(gca,'ytick',[],'xtick',[0 cumsum(SegDurs)'],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
box off
axis fill

suptitle('Stim concatenated')

% Save
print_eps_kp(gcf,fullfile(savedir,'Cluster_ConcatStim'))


