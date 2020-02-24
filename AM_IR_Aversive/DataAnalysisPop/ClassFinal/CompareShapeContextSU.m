% Compare Shape vs Context  :  SU 

whichStim = 'Speech';
whichClass = 'OnlyTemp';

% Data settings
fn = set_paths_directories('','',1);

switch whichStim
    case 'AC'    
        keyboard
        
        q = load(fullfile(fn.processed,'Units'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % SHAPE: Load SU results
        datadir = fullfile(fn.figs,'ClassAM','AC',whichClass,'each');
        q=load(fullfile(datadir,'CR_each.mat'));
        CRe_S_Full = q.CR;
        clear q
        
        % CONTEXT: load SU results
        datadir = fullfile(fn.figs,'ClassContext','AC',whichClass,'each');
        q=load(fullfile(datadir,'CR_each.mat'));
        CRe_Context = q.CR;
        clear q
        
        % Also get CTTS
        [~,theseCells] = recallDataParams('AC','each');
        
        savedir = fullfile(fn.figs,'ClassContext','AC','CompareShape');
        
    case 'Speech'        
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitInfo = q.UnitInfo;
        clear q
        
        % SHAPE: Load SU results
        
        % Full
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','Full','each','CR_each.mat'));
        CRe_S_Full = q.CR;
        clear q
        % Temp
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyTemp','each','CR_each.mat'));
        CRe_S_Temp = q.CR;
        clear q
        % Rate
        q=load(fullfile(fn.figs,'ClassSpeech','Speech','OnlyRate','each','CR_each.mat'));
        CRe_S_Rate = q.CR;
        clear q
        
        % CONTEXT: load SU results
        
        % Full
        q=load(fullfile(fn.figs,'ClassContext','Speech','Full','each','CR_each.mat'));
        CRe_C_Full = q.CR;
        clear q
        % Temp
        q=load(fullfile(fn.figs,'ClassContext','Speech','OnlyTemp','each','CR_each.mat'));
        CRe_C_Temp = q.CR;
        clear q
        % Rate
        q=load(fullfile(fn.figs,'ClassContext','Speech','OnlyRate','each','CR_each.mat'));
        CRe_C_Rate = q.CR;
        clear q
        
        % Also get CTTS
        [~,theseCells] = recallDataParams('Speech','each');
        
        savedir = fullfile(fn.figs,'ClassContext','Speech','CompareShape');
end
if ~exist(savedir,'dir')
    mkdir(savedir)
end


% Fig settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
% fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallsmall = [1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2];
widesmall = [1 scrsz(4)/3 scrsz(3)/3*2 scrsz(4)/3];




%% For each Class Type, compare Shape and Context SU distributions

ClassTypes = {'Full' 'OnlyTemp' 'OnlyRate'};

for iwc = 1:numel(ClassTypes)
    
    whichClass = ClassTypes{iwc};
    switch whichClass
        case 'Full'
            CR_S = CRe_S_Full;
            CR_C = CRe_C_Full;
        case 'OnlyTemp'
            CR_S = CRe_S_Temp;
            CR_C = CRe_C_Temp;
        case 'OnlyRate'
            CR_S = CRe_S_Rate;
            CR_C = CRe_C_Rate;
    end
    
    
    dp_edges = -0.5:0.1:3.5;
    dp_xvals = dp_edges(1:end-1)+mode(diff(dp_edges))/2;
    
    % Shape first
    % plot and measure mean & skewness
    [dpShape,iCReS] = sort(CR_S.dprime,'descend');
    S = skewness(CR_S.dprime,0);
    
    H_Shape = histcounts(dpShape,dp_edges,'Normalization','probability');
    
    
    % Context next
    % one seg at a time
    Segs = unique(CR_C.Seg)';
    
    dpContextSeg = nan(size(UnitInfo,1),numel(Segs));
    iCReCSeg     = nan(size(UnitInfo,1),numel(Segs));
    S_Seg        = nan(1,numel(Segs));
    
    H_Context    = nan(numel(dp_xvals),numel(Segs));
    for is = Segs
        [dpContextSeg(1:sum(CR_C.Seg==is),is),iCReCSeg(1:sum(CR_C.Seg==is),is)]...
            = sort(CR_C.dprime(CR_C.Seg==is),'descend');
        S_Seg(1,is) = skewness(dpContextSeg(1:sum(CR_C.Seg==is),is),0);
        
        H_Context(:,is) = histcounts(CR_C.dprime(CR_C.Seg==is),dp_edges,...
            'Normalization','probability');
        
    end
    
    
    figure;
    set(gcf,'Position',widesmall)
    plot(dp_xvals,H_Shape,'k','LineWidth',2)
    hold on
    plot(dp_xvals,mean(H_Context,2),'b','LineWidth',2)
    
    savename = sprintf('dpPDF_ShapeContext_%s',whichClass);
    print_eps_kp(gcf,fullfile(savedir,savename))
    
end



%%  Removing temporal info has less effect on Context than Shape


figure;
set(gcf,'Position',widesmall)

subplot(1,2,1)
plot([-0.5 3],[-0.5 3],'-','Color',0.7*[1 1 1])
hold on
plot(CRe_S_Full.dprime,CRe_S_Rate.dprime,'.k')
plot(CRe_C_Full.dprime,CRe_C_Rate.dprime,'.m')
plot(mean(CRe_S_Full.dprime),mean(CRe_S_Rate.dprime),'o','MarkerSize',10,'LineWidth',3,'MarkerFaceColor',[1 0.99 1],'MarkerEdgeColor',[0 0 0.01])
plot(mean(CRe_C_Full.dprime),mean(CRe_C_Rate.dprime),'o','MarkerSize',10,'LineWidth',3,'MarkerFaceColor','y','MarkerEdgeColor',[0 0 0.01])
axis square
xlim([-0.5 3])
ylim([-0.5 3])
ylabel('Rate only')
xlabel('Full classifier')
title(sprintf('Rate d''/Full d''\nShape: %0.2f, Context: %0.2f', mean(CRe_S_Rate.dprime)/mean(CRe_S_Full.dprime), mean(CRe_C_Rate.dprime)/mean(CRe_C_Full.dprime) ))

% mean(CRe_S_Rate.dprime)/mean(CRe_S_Full.dprime)
% mean(CRe_C_Rate.dprime)/mean(CRe_C_Full.dprime)


subplot(1,2,2)
plot([-0.5 3],[-0.5 3],'-','Color',0.7*[1 1 1])
hold on
plot(CRe_S_Full.dprime,CRe_S_Temp.dprime,'.k')
plot(CRe_C_Full.dprime,CRe_C_Temp.dprime,'.b')
plot(mean(CRe_S_Full.dprime),mean(CRe_S_Temp.dprime),'o','MarkerSize',10,'LineWidth',3,'MarkerFaceColor',[1 0.99 1],'MarkerEdgeColor',[0 0 0.01])
plot(mean(CRe_C_Full.dprime),mean(CRe_C_Temp.dprime),'o','MarkerSize',10,'LineWidth',3,'MarkerFaceColor','y','MarkerEdgeColor',[0 0 0.01])
axis square
xlim([-0.5 3])
ylim([-0.5 3])
ylabel('Temporal only')
xlabel('Full classifier')
title(sprintf('Temp d''/Full d''\nShape: %0.2f, Context: %0.2f', mean(CRe_S_Temp.dprime)/mean(CRe_S_Full.dprime), mean(CRe_C_Temp.dprime)/mean(CRe_C_Full.dprime) ))

% mean(CRe_S_Temp.dprime)/mean(CRe_S_Full.dprime)
% mean(CRe_C_Temp.dprime)/mean(CRe_C_Full.dprime)

savename = sprintf('dp_FullRateTemp_SU');
print_eps_kp(gcf,fullfile(savedir,savename))




keyboard








