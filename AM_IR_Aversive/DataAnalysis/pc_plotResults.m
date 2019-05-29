function pc_plotResults
% 

%% 

% Load Unit files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
[sigUnits,UnitData] = identifyResponsiveUnits(UnitData);

% Load classifier results
load(fullfile(fn.processed,'MPHclassifier','ClassData'));
q=load(fullfile(fn.processed,'MPHclassifier','ClassData_shuff'));
DataSh = q.Data;
clear q

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
scrsz = get(0,'ScreenSize');  %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
longrect   = [1 scrsz(4) scrsz(3) scrsz(4)/4];
halfrect   = [1 scrsz(4) scrsz(3) scrsz(4)/2];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [  84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;

% Other settings
histbinsize = 0.025;
AMrates = [2 4 8 16 32];

% Plot
hf=figure; 
set(hf,'Position',fullscreen,'NextPlot','add')
hold on


%%  OVERALL DISTRIBUTION OF D' VALUES 

alfa = 0.05;

N  = 0;
Ns = 0;
Nb = 0;
plot_vals  = [];
BMF_vals   = [];
plot_v_Sh  = [];
BMF_v_Sh   = [];
sig_vals   = [];
nsr_vals   = [];

for iUn = 1:size(Data,1)
    
    this_BMF = UnitData(iUn).iBMF_FR;
    
    for irate = 1:size(Data,2)
        
        N = N + size(Data(iUn,irate).Res_L1o.dprime,1);
        
        plot_vals = [plot_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
        plot_v_Sh = [plot_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
        
        if UnitData(iUn).iBMF_FR == irate
            BMF_vals = [BMF_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
            BMF_v_Sh = [BMF_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
            Nb = Nb + size(Data(iUn,irate).Res_L1o.dprime,1);
        end
        
        if UnitData(iUn).kw_p<alfa && UnitData(iUn).wx_p<alfa
            sig_vals = [sig_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
            Ns = Ns + size(Data(iUn,irate).Res_L1o.dprime,1);
        else 
            nsr_vals = [nsr_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
        end
    end
    
end

plot_vals(plot_vals<0) = 0;
BMF_vals(BMF_vals<0)   = 0;
plot_v_Sh(plot_v_Sh<0) = 0;
BMF_v_Sh(BMF_v_Sh<0)   = 0;
sig_vals(sig_vals<0)   = 0;
nsr_vals(nsr_vals<0)   = 0;


% Plot
subplot(2,2,1); hold on
plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
plot([2 2],[0 1],'--','Color',0.5.*[1 1 1])
plot([3 3],[0 1],'--','Color',0.5.*[1 1 1])
ih(1)=histogram(plot_v_Sh,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','cdf');
ih(2)=histogram(plot_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','k','LineWidth',4,'Normalization','cdf');
ih(3)=histogram(sig_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','b','LineWidth',2,'Normalization','cdf');
ih(4)=histogram(nsr_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','c','LineWidth',2,'Normalization','cdf');
ih(5)=histogram(BMF_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',[181   0  52]./255,'LineWidth',2,'Normalization','cdf');
xlabel('d prime')
ylabel('Cumulative probability')
ylim([0 1])
xlim([0 2])
legend(ih,{'shuff'  sprintf('ALL (%i)',N)   sprintf('sigTC (%i)',Ns)  'nsTC'  sprintf('BMF (%i)',Nb) },'Location','southeast')
title('d'' distribution for all MPH comparisons')



%%  HISTOGRAM OF D' VALUES BY AM RATE

clear ih
subplot(2,2,2); hold on

for irate = 1:size(Data,2)
    
    dp_plot = [];
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data)
            continue
        end
        
        dp_plot  = [dp_plot Data(iUn,irate).Res_L1o.dprime(:,2)'];
        dp_plot(dp_plot<0) = 0;
    end
    
    ih(irate)=histogram(dp_plot,0:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor',colors(irate,:),'LineWidth',2,'Normalization','cdf');
    
end
legend(ih,cellfun(@num2str,num2cell(AMrates),'UniformOutput',false),'Location','southeast')
xlim([0 2])
ylim([0 1])
xlabel('d prime')

title('d'' by AM rate')



%%  HISTOGRAM OF D' VALUES BY RELATION TO BMF

Nb=0;
plot_vals = [];
clear ih

for iUn = 1:size(Data,1)
    if isempty(UnitData(iUn).iBMF_FR)
        continue
    else 
        Nb = Nb+1;
    end
    
    for irate = 1:size(Data,2)
        dBMF      = irate - UnitData(iUn).iBMF_FR;
        dp_plot  =  Data(iUn,irate).Res_L1o.dprime(:,2);
        dp_plot(dp_plot<0) = 0;
        plot_vals = [plot_vals; dBMF.*ones(size(dp_plot)) dp_plot];
    end
    
end

% Plot
subplot(2,2,4); hold on

legstr = cell(1,9);
colors_bmf = cmocean('delta',13);
colors_bmf([1 2 12 13],:) = [];
colors_bmf(5,:) = [0 0 0];

for isp = 1:9
    
    iBMF = plot_vals(:,1) == (isp - 5);
    
    ih(isp)=histogram(plot_vals(iBMF,2),0:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor',colors_bmf(isp,:),'LineWidth',2,'Normalization','cdf');
    
    legstr{isp} = ['BMF + ' num2str(isp - 5)];
end

legend(ih,legstr,'Location','southeast')

xlim([0 2])
ylim([0 1])
xlabel('d prime')

title(sprintf('d'' by this period''s log distance to BMF (%i units with BMF)',Nb))


suptitle(sprintf('MPH Context Classifier Results\nN=%i units',size(Data,1)))


% Save 
print_eps_kp(hf,fullfile(fn.figs,'MPHclass','dprimeDistributions'))


keyboard


%%  Compare performance of classifiers

dp_L1o = [];
dp_t10 = [];

for iUn = 1:size(Data,1)
    
    for irate = 1:size(Data,2)
        dp_L1o  =  [dp_L1o; Data(iUn,irate).Res_L1o.dprime(:,2)];
        dp_t10  =  [dp_t10; Data(iUn,irate).Res_t10.dprime(:,2)];
    end
    
end
dp_L1o(dp_L1o<0) = 0;
dp_t10(dp_t10<0) = 0;


% Plot
hfcc = figure;
plot([0 4]',[0 4]','-k')
hold on
plot(dp_L1o,dp_t10,'ok','LineWidth',1)

axis square
xlabel('dprime: Leave 1 out')
ylabel('dprime: 10 trial templates')

title('Comparing classifier results')

% Save 
print_eps_kp(hfcc,fullfile(fn.figs,'MPHclass','CompareClassifiers'))




end


