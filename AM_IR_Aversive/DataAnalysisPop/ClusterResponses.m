



%%%%  mark or separate NS cells 



% close all
useFR     =   'log'; 
whichStim =   'Speech';

savename = sprintf('Cluster_%s_%s_ConcatStim',whichStim,useFR);


% Load Unit data files
fn = set_paths_directories('','',1);
switch whichStim
    
    case 'AM'
        
        q = load(fullfile(fn.processed,'Units'));
        
        % Load MPH data
        load(fullfile(fn.figs,'PopMPH','MPHdata_gauss'));
        
        % Filter to just 2, 4, 8 Hz
        
        
        SegDurs  = [500; 250; 125; 63; 32];

    
    case 'Speech'
        
        q = load(fullfile(fn.processed,'UnitsVS'));
        
        % Load RepSpeech templates
        k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
        kfns = fieldnames(k);
        
        Duration = max(structfun(@length,k));
        Segs     = 1:numel(structfun(@length,k));
        SegDurs  = structfun(@length,k);
        
        load(fullfile(fn.figs,'PopulationSpeechSegments','RepSpeech'))

end
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end

savedir = fullfile(fn.figs,'ClusterResponses');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Fig settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',8)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


%% Data settings

switch useFR
    case 'z'
        ThisData    = zFR_vec;
    case 'log'
        ThisData    = FR_vec;
end



%% Main 

iUn_GoodData = find(sum(sum(~isnan(ThisData),2),3)==sum(SegDurs));

GoodData = [];
for ist = 1:size(ThisData,3)
    
    pdms = SegDurs(ist);
    
    add_data = ThisData(iUn_GoodData,1:pdms,ist);
%     add_data = repmat(add_data,[1 ceil(500/pdms)]);
    
    GoodData = [GoodData add_data];
    
end
% clipVal  = 10000;
% GoodData(GoodData> clipVal)  =  clipVal;
% GoodData(GoodData<-clipVal)  = -clipVal;


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

print_eps_kp(gcf,fullfile(savedir,savename))


