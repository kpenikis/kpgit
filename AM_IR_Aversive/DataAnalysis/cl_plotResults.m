function cl_plotResults
%  
%

%%%%%%%%%%%%%%%%%%
FR_cutoff = 0;
%%%%%%%%%%%%%%%%%%
dp_cutoff = 0;
%%%%%%%%%%%%%%%%%%
alfa      = 0.05; 
%%%%%%%%%%%%%%%%%%
PLOT_F1   = 1;
%%%%%%%%%%%%%%%%%%

% Load Unit files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;     clear q
[sigUnits,UnitData] = identifyResponsiveUnits(UnitData);

% Load classifier results
load(fullfile(fn.processed,'MPHclassifier','ClassData'));
q=load(fullfile(fn.processed,'MPHclassifier','ClassData_shuff'));
DataSh = q.Data;           clear q

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
scrsz = get(0,'ScreenSize');  %[left bottom width height]
halfscreen = [1 scrsz(4) scrsz(3) scrsz(4)/2];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [  84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
% Other settings
histbinsize = 0.0001;
AMrates = [2 4 8 16 32];


if PLOT_F1
%%  FIGURE 1: HISTOGRAMS

% Create plot
hf1=figure; 
set(hf1,'Position',halfscreen,'NextPlot','add')
hold on


%---------------------------------------
%  OVERALL DISTRIBUTION OF D' VALUES 
%---------------------------------------

N  = 0;    Ns = 0;   Nb = 0;
plot_vals  = [];
BMF_vals   = [];
plot_v_Sh  = [];
BMF_v_Sh   = [];
sig_vals   = [];
nsr_vals   = [];
dp_facil   = [];
dp_adapt   = [];
dp_nchng   = [];
PdcRespType = cell(size(DataSh,1),5);


for iUn = 1:size(DataSh,1)
    
    this_BMF = UnitData(iUn).iBMF_FR;
    
    for irate = 1:size(Data,2)
        
        if isempty(Data(iUn,irate).data)
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
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
        
        
        % Get Pdc response type                      [ early  late  p_val ]
        PdcRespType{iUn,irate} = 'na';
        if ~isempty(UnitData(iUn).DeltaNspk) && ~isempty(UnitData(iUn).DeltaNspk{irate})
            if (UnitData(iUn).DeltaNspk{irate}(3))<alfa && (UnitData(iUn).DeltaNspk{irate}(1) > UnitData(iUn).DeltaNspk{irate}(2))
                % adapting
                PdcRespType{iUn,irate} = 'A';
                dp_adapt = [dp_adapt; Data(iUn,irate).Res_L1o.dprime(:,2)];
            elseif (UnitData(iUn).DeltaNspk{irate}(3))<alfa && (UnitData(iUn).DeltaNspk{irate}(1) < UnitData(iUn).DeltaNspk{irate}(2))
                % facilitating
                PdcRespType{iUn,irate} = 'F';
                dp_facil = [dp_facil; Data(iUn,irate).Res_L1o.dprime(:,2)];
            else
                PdcRespType{iUn,irate} = 'NC';
                dp_nchng = [dp_nchng; Data(iUn,irate).Res_L1o.dprime(:,2)];
            end
        end
        
        
        
    end %irate
    
end %iUn

plot_vals(plot_vals<0) = 0;
BMF_vals(BMF_vals<0)   = 0;
plot_v_Sh(plot_v_Sh<0) = 0;
BMF_v_Sh(BMF_v_Sh<0)   = 0;
sig_vals(sig_vals<0)   = 0;
nsr_vals(nsr_vals<0)   = 0;
dp_adapt(dp_adapt<0)   = 0;
dp_facil(dp_facil<0)   = 0;
dp_nchng(dp_nchng<0)   = 0;


%---------------------------------------
% Overall distributions of d'
%---------------------------------------
subplot(1,3,1); hold on
grid on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
% plot([2 2],[0 1],'--','Color',0.5.*[1 1 1])
% plot([3 3],[0 1],'--','Color',0.5.*[1 1 1])
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
axis square
legend(ih,{'shuff'  sprintf('ALL (%i)',N)   sprintf('sigTC (%i)',Ns)  'nsTC'  sprintf('BMF (%i)',Nb) },'Location','southeast')
title('d'' distribution for all MPH comparisons')

% Stats
h1 = adtest(plot_vals);
h2 = adtest(plot_v_Sh);

p_rs = ranksum(plot_vals,plot_v_Sh);
[h_ks,p_ks] = kstest2(plot_vals,plot_v_Sh);


%---------------------------------------
%  Distribution of d' by Pdc resp type
%---------------------------------------
clear ih

subplot(1,3,2); hold on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
grid on
ih(1)=histogram(dp_nchng,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.1.*[1 1 1],'LineWidth',2,'Normalization','cdf');
ih(2)=histogram(dp_adapt,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[1 1 0],'LineWidth',2,'Normalization','cdf');
ih(3)=histogram(dp_facil,0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[0 1 0],'LineWidth',2,'Normalization','cdf');
ih(4)=histogram([dp_facil; dp_adapt],0:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[0 0 1],'LineWidth',2,'Normalization','cdf');

legend(ih,{['NC (n=' num2str(length(dp_nchng)) ')'] ['Adapting (n=' num2str(length(dp_adapt)) ')'] ['Facilitating (n=' num2str(length(dp_facil)) ')'] ['A or F (n=' num2str(length(dp_adapt)+length(dp_facil)) ')']},'Location','southeast')
xlim([0 2])
ylim([0 1])
xlabel('d prime')
axis square
title('d'' by Pdc response type')

% Stats

[h_ks_2,p_ks_2] = kstest2(dp_nchng(round(linspace(1,length(dp_nchng),200))),[dp_facil; dp_adapt]);
% [h_ks_2,p_ks_2] = kstest2(dp_nchng,[dp_facil; dp_adapt]);


%---------------------------------------
%  Distribution of d' by AM rate
%---------------------------------------
clear ih
legstr = cell(1,size(Data,2));

subplot(1,3,3); hold on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
grid on

for irate = 1:size(Data,2)
    
    dp_plot = [];
    N = 0;
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data)
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
        dps = Data(iUn,irate).Res_L1o.dprime(:,2)';
        dps(isnan(dps)) = [];
        dps(dps<0) = 0;
        
        dp_plot  = [dp_plot dps];
        N        = N+size(Data(iUn,irate).Res_L1o.dprime,1); %length(dps);
    end
    
    ih(irate)=histogram(dp_plot,0:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor',colors(irate,:),'LineWidth',2,'Normalization','cdf');
    
    legstr{irate} = sprintf('%i Hz (%i)',AMrates(irate),N);

end

legend(ih,legstr,'Location','southeast')
xlim([0 2])
ylim([0 1])
xlabel('d prime')
axis square

title('d'' by AM rate')

suptitle(sprintf('MPH Context Classifier Results\nN=%i units',size(Data,1)))


% Save 
print_eps_kp(hf1,fullfile(fn.figs,'MPHclass','dprimeDistributions'))


%%  
end

%%  FIGURE 2: CORRELATIONS

hf2=figure;
set(hf2,'Position',fullscreen,'NextPlot','add')
hold on


%%  Relationship of d' values to BMF 

Nb=0;
plot_vals = [];
clear ih

for iUn = 1:size(Data,1)
    
    if isempty(UnitData(iUn).iBMF_FR)
        continue
    end
    
    Nb = Nb+1;
    
    for irate = 1:size(Data,2)
        
        if isempty(Data(iUn,irate).data)
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
        dBMF      = irate - UnitData(iUn).iBMF_FR;
        dp_plot  =  Data(iUn,irate).Res_L1o.dprime(:,2);
        dp_plot(dp_plot<0) = 0;
        plot_vals = [plot_vals; dBMF.*ones(size(dp_plot)) dp_plot];
    end
    
end


% Plot  (beeswarm)
idx = plot_vals(:,2)>1;

subplot(2,3,1); hold on
plotSpread( plot_vals(:,2), 'distributionIdx', plot_vals(:,1), 'distributionColors', 'k', 'showMM', 3 );

xlabel('log2 distance to BMF')
ylabel('dprime')
title(sprintf('Relationship of d'' to AM rate tuning (%i units with BMF)',Nb))

% Plot (cdf)
% subplot(2,3,1); hold on
% 
% legstr = cell(1,9);
% colors_bmf = cmocean('delta',13);
% colors_bmf([1 2 12 13],:) = [];
% colors_bmf(5,:) = [0 0 0];
% 
% for isp = 1:9
%     iBMF = plot_vals(:,1) == (isp - 5);
%     ih(isp)=histogram(plot_vals(iBMF,2),0:histbinsize:4,'DisplayStyle','stairs',...
%         'EdgeColor',colors_bmf(isp,:),'LineWidth',2,'Normalization','cdf');
%     legstr{isp} = ['BMF + ' num2str(isp - 5)];
% end
% 
% legend(ih,legstr,'Location','southeast')
% xlim([0 2])
% ylim([0 1])
% xlabel('d prime')
% title(sprintf('d'' by this period''s log distance to BMF (%i units with BMF)',Nb))


%%  Relationship to trial variability 

plotvals_FR      = [];
plotvals_FF      = [];
plotvals_FF_best = [];
plotvals_CS      = [];
plotvals_CS_best = [];

% plotvals_FF_bmf  = [];
% plotvals_CS_bmf  = [];

Nb=0;
clear ih

for iUn = 1:size(Data,1)
    
%     if isempty(UnitData(iUn).iBMF_FR)
%         continue
%     end
    
    FR_pdc = nan(1,size(Data,2));
    FF_pdc = nan(1,size(Data,2));
    CS_pdc = nan(1,size(Data,2));
    
    for irate = 1:4
        if isempty(Data(iUn,irate).data)
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
        % d prime
        dp_plot  =  Data(iUn,irate).Res_L1o.dprime(:,2);
        dp_plot(dp_plot<dp_cutoff) = [];
        dp_plot(isnan(dp_plot))    = [];
        if isempty(dp_plot)
            continue
        end
        
        % FR
        FR_pdc(irate) = mean(mean(Data(iUn,irate).data(1).raster,2));
        plotvals_FR = [plotvals_FR; FR_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
        % Fano factor
        FF_pdc(irate) = var(mean(Data(iUn,irate).data(1).raster,2)) / mean(mean(Data(iUn,irate).data(1).raster,2));
        
        if any(FF_pdc<0.005)
            aaa=234;
        end
        
        plotvals_FF = [plotvals_FF; FF_pdc(irate).*ones(length(dp_plot),1) dp_plot];
%         plotvals_FF_bmf = [plotvals_FF_bmf; FF_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
        % Corr score 
        nt = min( size(Data(iUn,irate).data(1).raster,1), 100 );
        c = nan(nt);
        for ii = 1:nt  % for all trials
            for j = 1:ii-1  % up until that trial
                c(ii,j) = corr( Data(iUn,irate).data(1).raster(ii,:)' , Data(iUn,irate).data(1).raster(j,:)' ); 
            end
        end
        CS_pdc(irate) = nanmedian(c(~isnan(c)));
        
        plotvals_CS = [plotvals_CS; CS_pdc(irate).*ones(length(dp_plot),1) dp_plot];
%         plotvals_CS_bmf = [plotvals_CS_bmf; CS_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
    end
    
    if sum(isnan(FF_pdc))>4
        continue
    end
    
    % Now get values only for best AM rate, defined by lowest variance in 
    % periodic context
    ir = [];
    [FFr,ir] = min(FF_pdc);
    
    Nb = Nb+1;
    
    dp_plot  =  Data(iUn,ir).Res_L1o.dprime(:,2);
    dp_plot(dp_plot<dp_cutoff) = [];
    plotvals_FF_best = [plotvals_FF_best; FF_pdc(ir).*ones(sum(~isnan(dp_plot)),1) dp_plot(~isnan(dp_plot))];
    
    
    % Now get values only for best AM rate, defined by lowest variance in
    % periodic context
    ir = [];
    [CSr,ir] = max(CS_pdc);
        
    dp_plot  =  Data(iUn,ir).Res_L1o.dprime(:,2);
    dp_plot(dp_plot<dp_cutoff) = [];
    plotvals_CS_best = [plotvals_CS_best; CS_pdc(ir).*ones(sum(~isnan(dp_plot)),1) dp_plot(~isnan(dp_plot))];
    
end

% FF (sp2)
Plot_Val_str = 'plotvals_FF';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,2); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('FF')


% FF best (sp3)
Plot_Val_str = 'plotvals_FF_best';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,3); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('FF')


% CS (sp5)
Plot_Val_str = 'plotvals_CS';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,5); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('Corr Score')


% CS best (sp6)
Plot_Val_str = 'plotvals_CS_best';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,6); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('Corr Score')


% FR (sp4)
Plot_Val_str = 'plotvals_FR';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,4); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('FR (sp/s)')
ylabel('dprime')



%% Finish and save

suptitle(sprintf('What underlies classifier results?\nN=%i units, FR cutoff %i, excluding 32 Hz MPH',size(Data,1),FR_cutoff))

print_eps_kp(gcf,fullfile(fn.figs,'MPHclass','dprimeRelationships_NewClassData_1'))


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


