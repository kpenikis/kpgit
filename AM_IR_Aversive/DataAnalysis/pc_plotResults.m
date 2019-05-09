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


% Other settings
histbinsize = 0.025;
AMrates = [2 4 8 16 32];



%%  OVERALL DISTRIBUTION OF D' VALUES 

N = 0;
plot_vals  = [];
BMF_vals   = [];
plot_v_Sh  = [];
BMF_v_Sh   = [];

for iUn = 1:size(Data,1)
    
    N = N+1;
    this_BMF = UnitData(iUn).iBMF_FR;
    
    for irate = 1:size(Data,2)
        
        plot_vals = [plot_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
        plot_v_Sh = [plot_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
        
        if UnitData(iUn).iBMF_FR == irate
            BMF_vals = [BMF_vals; Data(iUn,irate).Res_L1o.dprime(:,2)];
            BMF_v_Sh = [BMF_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
        end
    end
    
end

plot_vals(plot_vals<0) = 0;
BMF_vals(BMF_vals<0)   = 0;
plot_v_Sh(plot_v_Sh<0) = 0;
BMF_v_Sh(BMF_v_Sh<0)   = 0;


% Plot
hf=figure; 
set(hf,'Position',halfrect,'NextPlot','add')
hold on

subplot(1,2,1); hold on
ih(1)=histogram(plot_v_Sh,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','probability');
ih(2)=histogram(plot_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','k','LineWidth',2,'Normalization','probability');
ih(3)=histogram(BMF_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',[181   0  52]./255,'LineWidth',2,'Normalization','probability');
legend(ih,{'shuff' 'all' 'just BMF'})
xlabel('d prime')
ylabel('Probability')


subplot(1,2,2); hold on
plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
plot([2 2],[0 1],'--','Color',0.5.*[1 1 1])
plot([3 3],[0 1],'--','Color',0.5.*[1 1 1])
ih(1)=histogram(plot_v_Sh,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','cdf');
ih(2)=histogram(plot_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','k','LineWidth',2,'Normalization','cdf');
ih(3)=histogram(BMF_vals,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',[181   0  52]./255,'LineWidth',2,'Normalization','cdf');
legend(ih,{'shuff' 'all' 'just BMF'})
xlabel('d prime')
ylabel('Cumulative probability')
ylim([0 1])

suptitle(sprintf('Discriminability of MPH from Periodic\nN=%i',N))

% Save 
print_eps_kp(hf,fullfile(fn.figs,'MPHclass','dprimeDistr'))




%%  HISTOGRAM OF D' VALUES BY AM RATE

hf=figure; 
set(hf,'Position',tallrect,'NextPlot','add')
hold on

N = 0;
for irate = 1:size(Data,2)
    
    dp_plot = [];
    dp_plot_s  = [];
    dp_plot_ns = [];
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data)
            continue
        end
        
        N = N+1;
        
        if UnitData(iUn).kw_p < 0.05
            dp_plot_s  = [dp_plot_s  Data(iUn,irate).Res_L1o.dprime(:,2)'];
            dp_plot_s(dp_plot_s<0) = 0;
        else 
            dp_plot_ns = [dp_plot_ns Data(iUn,irate).Res_L1o.dprime(:,2)'];
            dp_plot_ns(dp_plot_ns<0) = 0;
        end
        dp_plot  = [dp_plot Data(iUn,irate).Res_L1o.dprime(:,2)'];
        dp_plot(dp_plot<0) = 0;
    end
    
    isp(irate)=subplot(size(Data,2),1,irate); hold on
    ih(1)=histogram(dp_plot_ns,0:0.25:4,'DisplayStyle','stairs',...
        'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','probability');
    ih(2)=histogram(dp_plot_s,0:0.25:4,'DisplayStyle','stairs',...
        'EdgeColor','k','LineWidth',2,'Normalization','probability');
    title([num2str(AMrates(irate)) ' Hz'])
end

linkaxes(isp,'xy')
xlim([0 4])
ylim([0 0.65])
xlabel('d prime')
ylabel('Prob')
legend(ih,{'non-sig resp' 'Sig resp'})


suptitle(sprintf('Discriminability of MPH from Periodic\nN=%i',N/size(Data,2)))

% Save 
print_eps_kp(hf,fullfile(fn.figs,'MPHclass','dprimeDistr_AMrate'))




%%  HISTOGRAM OF D' VALUES BY RELATION TO BMF

N = 0;
plot_vals = [];

for iUn = 1:size(Data,1)
    if isempty(UnitData(iUn).iBMF_FR)
        continue
    else 
        N = N+1;
    end
    
    for irate = 1:size(Data,2)
        dBMF      = irate - UnitData(iUn).iBMF_FR;
        dp_plot  =  Data(iUn,irate).Res_L1o.dprime(:,2);
        dp_plot(dp_plot<0) = 0;
        plot_vals = [plot_vals; dBMF.*ones(size(dp_plot)) dp_plot];
    end
    
end

% Plot
hf=figure; 
set(hf,'Position',longrect,'NextPlot','add')
hold on
for isp = 1:9
    
    iBMF = plot_vals(:,1) == (isp - 5);
    
    hsp(isp)=subplot(1,9,isp); hold on
    
    histogram(plot_vals(iBMF,2),0:0.25:4,'DisplayStyle','stairs',...
        'EdgeColor','k','LineWidth',2,'Normalization','probability');
    plot(median(plot_vals(iBMF,2)),0,'^g','LineWidth',2)
    
    title(['BMF + ' num2str(isp - 5)])
end
linkaxes(hsp,'xy')
xlim([0 4])
ylim([0 0.4])
xlabel('d prime')
ylabel('Prob')

suptitle(sprintf('Discriminability of MPH from Periodic  |  N=%i with BMF',N))

% Save 
print_eps_kp(hf,fullfile(fn.figs,'MPHclass','dprimeDistr_dBMF'))





%%  Compare performance of classifiers

N = 0;
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





%%  FR HISTORY DEPENDENCE FOR HIGH D' 

dp_thresh = 1.5;
thisMet = 'Prev100msFR';

x1_plot  = [];
x2_plot  = [];

hf=figure; 
set(hf,'Position',fullscreen,'NextPlot','add')

Nh = 0;
Nl = 0;
for irate = 1:size(Data,2)
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data), continue, end
        
        dprimes = Data(iUn,irate).Res_L1o.dprime(:,2);
        high_dp = find(dprimes>dp_thresh)';
        low_dp  = find(dprimes<=dp_thresh)';
        
        for idp = low_dp
            
            Nl = Nl+1;
            
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);            
            
            x1_plot   = Data(iUn,irate).data(1).(thisMet);
            x2_plot   = Data(iUn,irate).data(id).(thisMet);
            
            figure(hf); hold on
            plot(x1_plot,x2_plot,'o','Color',0.7*[1 1 1],'LineWidth',2)
%             plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'o','Color',0.5*[1 1 1],'LineWidth',2)
            
        end %low_dp
        for idp = high_dp
            
            Nh = Nh+1;
            
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);            
            
            x1_plot   = Data(iUn,irate).data(1).(thisMet);
            x2_plot   = Data(iUn,irate).data(id).(thisMet);
            
            figure(hf); hold on
            plot(x1_plot,x2_plot,'ob','LineWidth',2)
            plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'ok','LineWidth',2)
            
        end %high_dp
        
    end %iUn
end %irate


% Finish and save fig 1
figure(hf); hold on
title(['Change in response as a function of ' thisMet ', N=' num2str(Nh) ', dp thresh=' num2str(dp_thresh)])
set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
axis square
plot([10^-2 10^2],[10^-2 10^2],'-c')

% print_eps_kp(hf,fullfile(fn.figs,'MPHclass',['RespDependence_dpAbove_' num2str(dp_thresh*10) '_' thisMet '_PreToPre']))



% Stat test for what datapoints are off unity line
% Take these dps, look at trial by trial FR history relationship 



%%  FR HISTORY DEPENDENCE FOR HIGH D' -- comparing change in prev/resp 

dp_thresh = 1.5;
thisMet = 'Prev100msFR';

x1_plot  = [];
x2_plot  = [];
y1_plot  = [];
y2_plot  = [];

hf=figure; 

hf2=figure;
set(hf2,'Position',longrect,'NextPlot','add')

N = 0;
for irate = 1:size(Data,2)
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data), continue, end
        
        dprimes = Data(iUn,irate).Res_L1o.dprime(:,2);
        high_dp = find(dprimes>dp_thresh)';
        
        for idp = high_dp
            
            N = N+1;
            
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);            
            dBMF      = irate - UnitData(iUn).iBMF_FR;
            
            x1_plot   = Data(iUn,irate).data(1).(thisMet);
            x2_plot   = Data(iUn,irate).data(id).(thisMet);
            
            y1_plot   = mean(mean(Data(iUn,irate).data(1).raster,1));
            y2_plot   = mean(mean(Data(iUn,irate).data(id).raster,1));
            
            figure(hf); hold on
%             plot([x1_plot x2_plot],[y1_plot y2_plot],'.-k')
            plot(x1_plot,x2_plot,'ok','LineWidth',2)
            plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'o','Color',0.5*[1 1 1],'LineWidth',2)
            
            if ~isempty(dBMF)
%                 figure(hf2); hold on
%                 hsp(dBMF+5)=subplot(1,9,dBMF+5); hold on
%                 plot([x1_plot x2_plot],[y1_plot y2_plot],'.-k')
%                 title(['BMF + ' num2str(dBMF)])
%                 set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
%                 axis square
            end
            
        end %high_dp
    end %iUn
end %irate


% Finish and save fig 1
figure(hf); hold on
title(['Change in response as a function of ' thisMet ', N=' num2str(N) ', dp thresh=' num2str(dp_thresh)])
set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
axis square
plot([10^-2 10^2],[10^-2 10^2],'-c')

print_eps_kp(hf,fullfile(fn.figs,'MPHclass',['RespDependence_dpAbove_' num2str(dp_thresh*10) '_' thisMet '_PreToPre']))


% Finish and save fig 2
% figure(hf2); hold on
% suptitle(['Change in response as a function of ' thisMet ', dp thresh=' num2str(dp_thresh)])
% 
% print_eps_kp(hf2,fullfile(fn.figs,'MPHclass',['RespDependence_dpAbove_' num2str(dp_thresh*10) '_' thisMet '_dBMF']))





end


