function cl_SU_plotResults
% cl_SU_plotResults
%   In the process of adapting for new data file structure, with matched MP
%   data separate from SU classifier results. Shuffled classifier results
%   have only 157 units still, while Res has 221. 
%

global AMrates

%%%%%%%%%%%%%%%%%%
FR_cutoff = 0;
%%%%%%%%%%%%%%%%%%
dp_cutoff = -inf;
%%%%%%%%%%%%%%%%%%
alfa      = 0.05; 
%%%%%%%%%%%%%%%%%%
PLOT_F1   = 0;
PLOT_F2   = 0;
PLOT_F3   = 1;
%%%%%%%%%%%%%%%%%%

% Load Unit files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;    clear q
[sigUnits,UnitData] = identifyResponsiveUnits(UnitData);
% Load MP data
q=load(fullfile(fn.processed,'MatchedMPs','MatchedMPdata'));
Data   = q.Data;          clear q
% Load classifier results
q=load(fullfile(fn.processed,'MPHclassifier','MPcontextSUclassRes'));
Res    = q.Res;           clear q
q=load(fullfile(fn.processed,'MPHclassifier','ClassData_shuff'));
DataSh = q.Data;          clear q


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
histbinsize = 0.01;


%%  FIGURE 1: HISTOGRAMS

if PLOT_F1

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
post_BMF   = [];
CellType   = {};
dp_broad   = [];
dp_narrow  = [];

bestResponders=find(cellfun(@max,{UnitData.FR_nrm})>quantile(cellfun(@max,{UnitData.FR_nrm}),0.85));
br_plotvals = [];

for iUn = 1:size(Data,1)
        
    this_BMF = UnitData(iUn).iBMF_FR;
    
    for irate = 1:size(Res,2)
        
        if isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'raster')
            continue
        elseif mean(mean(Data(iUn,irate).data(1,1).raster,2))<FR_cutoff
            continue
        end
        
        N = N + size(Res(iUn,irate).L1o.dprime,1);
        
        plot_vals = [plot_vals; Res(iUn,irate).L1o.dprime(:,2)];
%         plot_v_Sh = [plot_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
        
        if UnitData(iUn).FR_nrm(irate+1)>0.5 %ismember(iUn,bestResponders) 
            br_plotvals = [br_plotvals; Res(iUn,irate).L1o.dprime(:,2)];
        end
        
        if UnitData(iUn).iBMF_FR == irate
            BMF_vals = [BMF_vals; Res(iUn,irate).L1o.dprime(:,2)];
%             BMF_v_Sh = [BMF_v_Sh; DataSh(iUn,irate).Res_L1o.dprime(:,2)];
            Nb = Nb + size(Res(iUn,irate).L1o.dprime,1);
        end
        
        if UnitData(iUn).kw_p<alfa && UnitData(iUn).wx_p<alfa
            sig_vals = [sig_vals; Res(iUn,irate).L1o.dprime(:,2)];
            Ns = Ns + size(Res(iUn,irate).L1o.dprime,1);
        else 
            nsr_vals = [nsr_vals; Res(iUn,irate).L1o.dprime(:,2)];
        end
        
        if  ~isempty(UnitData(iUn).iBMF_FR) 
            
            data_idx = abs([Data(iUn,irate).data(:,2).PrevAMrt100] - AMrates(UnitData(iUn).iBMF_FR) * ones(size([Data(iUn,irate).data(:,2).PrevAMrt100]))) < 0.5;
            
            post_BMF = [post_BMF; Res(iUn,irate).L1o.dprime(data_idx,2)];
            
        end
        
        
        % Get Cell type           
        
        if UnitInfo(iUn,:).TroughPeak > 0.5
            % Broad
            CellType{end+1} = 'B';
            dp_broad  = [dp_broad; Res(iUn,irate).L1o.dprime(:,2)];
        else
            % Narrow
            CellType{end+1} = 'N';
            dp_narrow = [dp_narrow; Res(iUn,irate).L1o.dprime(:,2)];
        end
        
        
        % Get Pdc response type                      [ early  late  p_val ]
        
        if ~isempty(UnitData(iUn).DeltaNspk) && ~isempty(UnitData(iUn).DeltaNspk{irate})
            if (UnitData(iUn).DeltaNspk{irate}(3))<alfa && (UnitData(iUn).DeltaNspk{irate}(1) > UnitData(iUn).DeltaNspk{irate}(2))
                % adapting
                PdcRespType{iUn,irate} = 'A';
                dp_adapt = [dp_adapt; Res(iUn,irate).L1o.dprime(:,2)];
            elseif (UnitData(iUn).DeltaNspk{irate}(3))<alfa && (UnitData(iUn).DeltaNspk{irate}(1) < UnitData(iUn).DeltaNspk{irate}(2))
                % facilitating
                PdcRespType{iUn,irate} = 'F';
                dp_facil = [dp_facil; Res(iUn,irate).L1o.dprime(:,2)];
            else
                PdcRespType{iUn,irate} = 'S';
                dp_nchng = [dp_nchng; Res(iUn,irate).L1o.dprime(:,2)];
            end
        else
            PdcRespType{iUn,irate} = 'n';
            dp_nchng = [dp_nchng; Res(iUn,irate).L1o.dprime(:,2)];
        end
        
        
    end %irate
    
end %iUn


%% Create plot

hf1=figure; 
set(hf1,'Position',halfscreen,'NextPlot','add')
hold on

%---------------------------------------
%    Overall distributions of d'
%---------------------------------------
subplot(1,4,1); hold on
grid on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
% plot([2 2],[0 1],'--','Color',0.5.*[1 1 1])
% plot([3 3],[0 1],'--','Color',0.5.*[1 1 1])
ih(1)=histogram(plot_v_Sh,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.5.*[1 1 1],'LineWidth',2,'Normalization','cdf');
ih(2)=histogram(plot_vals,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','k','LineWidth',4,'Normalization','cdf');
ih(3)=histogram(sig_vals,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','b','LineWidth',2,'Normalization','cdf');
ih(4)=histogram(nsr_vals,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor','c','LineWidth',2,'Normalization','cdf');
ih(5)=histogram(BMF_vals,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',[181   0  52]./255,'LineWidth',2,'Normalization','cdf');
ih(6)=histogram(post_BMF,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',[240   140  0]./255,'LineWidth',2,'Normalization','cdf');
xlabel('d prime')
ylabel('Cumulative probability')
ylim([0 1])
xlim([-2 2])
axis square
legend(ih,{'shuff'  sprintf('ALL (%i)',N)   sprintf('sigTC (%i)',Ns)  'nsTC'  sprintf('BMF (%i)',Nb) sprintf('post-BMF (%i)',length(post_BMF)) },'Location','northwest')
title('d'' distribution for all MPH comparisons')

% Stats
h1 = adtest(plot_vals);
h2 = adtest(plot_v_Sh);

p_rs = ranksum(plot_vals,plot_v_Sh);
[h_ks,p_ks] = kstest2(plot_vals,plot_v_Sh);


%---------------------------------------
%  Distribution of d' by Cell type
%---------------------------------------
clear ih

subplot(1,4,2); hold on
grid on
ih(1)=histogram(plot_vals,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.1.*[1 1 1],'LineWidth',2,'Normalization','cdf');    %all
ih(2)=histogram(dp_broad,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.99.*[0.9 0.2 0.2],'LineWidth',2,'Normalization','cdf');    %broad
ih(3)=histogram(dp_narrow,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.99*[0.2 0.2 0.9],'LineWidth',2,'Normalization','cdf');    %narrow

legend(ih,{['All dps (n=' num2str(length(plot_vals)) ')'] ['Broad (n=' num2str(length(dp_broad)) ')'] ['Narrow (n=' num2str(length(dp_narrow)) ')'] },'Location','northwest')
xlim([-2 2])
ylim([0 1])
xlabel('d prime')
axis square
title('d'' by Cell type')


%---------------------------------------
%  Distribution of d' by Pdc resp type
%---------------------------------------
clear ih

subplot(1,4,3); hold on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
grid on
ih(1)=histogram(dp_nchng,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.1.*[1 1 1],'LineWidth',2,'Normalization','cdf');
ih(2)=histogram(dp_adapt,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[1 1 0],'LineWidth',2,'Normalization','cdf');
ih(3)=histogram(dp_facil,-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[0 1 0],'LineWidth',2,'Normalization','cdf');
ih(4)=histogram([dp_facil; dp_adapt],-4:histbinsize:4,'DisplayStyle','stairs',...
    'EdgeColor',0.9.*[0 0 1],'LineWidth',2,'Normalization','cdf');

legend(ih,{['NC (n=' num2str(length(dp_nchng)) ')'] ['Adapting (n=' num2str(length(dp_adapt)) ')'] ['Facilitating (n=' num2str(length(dp_facil)) ')'] ['A or F (n=' num2str(length(dp_adapt)+length(dp_facil)) ')']},'Location','northwest')
xlim([-2 2])
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

subplot(1,4,4); hold on
% plot([1 1],[0 1],'--','Color',0.5.*[1 1 1])
grid on

for irate = 1:size(Data,2)
    
    dp_plot = [];
    N = 0;
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'raster')
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
        dps = Res(iUn,irate).L1o.dprime(:,2)';
        dps(isnan(dps)) = [];
%         dps(dps<0) = 0;
        
        dp_plot  = [dp_plot dps];
        N        = N+size(Res(iUn,irate).L1o.dprime,1); %length(dps);
    end
    
    ih(irate)=histogram(dp_plot,-4:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor',colors(irate,:),'LineWidth',2,'Normalization','cdf');
    
    legstr{irate} = sprintf('%i Hz (%i)',AMrates(irate),N);

end

legend(ih,legstr,'Location','northwest')
xlim([-2 2])
ylim([0 1])
xlabel('d prime')
axis square

title('d'' by AM rate')

suptitle(sprintf('MPH Context Classifier Results\nN=%i units',size(Data,1)))


%% Save 
% print_eps_kp(hf1,fullfile(fn.figs,'MPHclass','dprimeDistributions'))

 
end %plot F1




%%  FIGURE 2: CORRELATIONS

if PLOT_F2
    
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
        dp_plot  =  Res(iUn,irate).L1o.dprime(:,2);
%         dp_plot(dp_plot<0) = 0;
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
plotvals_mad     = [];

% plotvals_FF_bmf  = [];
% plotvals_CS_bmf  = [];

Nb=0;
clear ih

for iUn = 1:size(Data,1)
    
%     if isempty(UnitData(iUn).iBMF_FR)
%         continue
%     end
    
    FR_pdc  = nan(1,size(Data,2));
    FF_pdc  = nan(1,size(Data,2));
    CS_pdc  = nan(1,size(Data,2));
    MAD_pdc = nan(1,size(Data,2));
    
    for irate = 1:4
        if isempty(Data(iUn,irate).data)
            continue
        elseif mean(mean(Data(iUn,irate).data(1).raster,2))<FR_cutoff
            continue
        end
        
        % d prime
        dp_plot  =  Res(iUn,irate).L1o.dprime(:,2);
        dp_plot(dp_plot<dp_cutoff) = [];
        dp_plot(isnan(dp_plot))    = [];
        if isempty(dp_plot)
            continue
        end
        
        %-------------
        %     FR
        FR_pdc(irate) = mean(mean(Data(iUn,irate).data(1).raster,2));
        plotvals_FR = [plotvals_FR; FR_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
        %-------------
        % Fano factor
        FF_pdc(irate) = var(mean(Data(iUn,irate).data(1).raster,2)) / mean(mean(Data(iUn,irate).data(1).raster,2));
        
        if any(FF_pdc<0.005)
            aaa=234;
        end
        
        plotvals_FF = [plotvals_FF; FF_pdc(irate).*ones(length(dp_plot),1) dp_plot];
%         plotvals_FF_bmf = [plotvals_FF_bmf; FF_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
        %-------------
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
        
        %-------------
        %     MAD
        MAD_pdc(irate) = mad(mean(Data(iUn,irate).data(1).raster,2),1);
%         MAD_irr(irate) = var(mean(Data(iUn,irate).data(1).raster,2)) / mean(mean(Data(iUn,irate).data(1).raster,2));
        
        plotvals_mad = [plotvals_mad; MAD_pdc(irate).*ones(length(dp_plot),1) dp_plot];
        
    end %irate
    
    if sum(isnan(FF_pdc))>4
        continue  % skip finding best 
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
%     ir = [];
%     [CSr,ir] = max(CS_pdc);
%         
%     dp_plot  =  Data(iUn,ir).Res_L1o.dprime(:,2);
%     dp_plot(dp_plot<dp_cutoff) = [];
%     plotvals_CS_best = [plotvals_CS_best; CS_pdc(ir).*ones(sum(~isnan(dp_plot)),1) dp_plot(~isnan(dp_plot))];
    
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
% Plot_Val_str = 'plotvals_CS';
Plot_Val_str = 'plotvals_mad';
Plot_Vals    = eval(Plot_Val_str);
Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];

subplot(2,3,5); hold on
axis square
plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')

[r,p]=corrcoef(Plot_Vals);
title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
xlabel('MAD')
% xlabel('Corr Score')


% CS best (sp6)
% Plot_Val_str = 'plotvals_CS_best';
% Plot_Vals    = eval(Plot_Val_str);
% Plot_Vals(sum(isnan(Plot_Vals),2)>0,:) = [];
% 
% subplot(2,3,6); hold on
% axis square
% plot(Plot_Vals(:,1),Plot_Vals(:,2),'ok')
% 
% [r,p]=corrcoef(Plot_Vals);
% title(sprintf( '%s  |  r=%0.2f, p=%0.1e', Plot_Val_str, r(1,2), p(1,2) )) 
% xlabel('Corr Score')


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

print_eps_kp(gcf,fullfile(fn.figs,'MPHclass','dprimeRelationships'))

end





%% FIGURE 3

% using VS as a proxy to identify "best responders"


if PLOT_F3

hf3=figure;
set(hf3,'Position',fullscreen,'NextPlot','add')
hold on

%---------------------------------------
%  Distribution of d' by AM rate
%---------------------------------------

for irate = 1:size(Data,2)
    
    dp_VS = []; 
    dp_nr = [];
    
    for iUn = 1:size(Data,1)
        if isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'raster')
            continue
        end
        
%         if (UnitData(iUn).VSdata_spk(2,irate+1) > 13.1) && (UnitData(iUn).VSdata_spk(1,irate+1) > 0.5)
        if (UnitData(iUn).VSdata_spk(2,irate+1) > (13.1 * 2.^(6-irate)) )
            dps = Res(iUn,irate).L1o.dprime(:,2)';
            dps(isnan(dps)) = [];            
            dp_VS  = [dp_VS dps];
            
        else
            dps = Res(iUn,irate).L1o.dprime(:,2)';
            dps(isnan(dps)) = [];
            dp_nr  = [dp_nr dps];
        end
    end %iUn
    
        
    subplot(2,3,irate); hold on
    grid on
    
    histogram(dp_nr,-4:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor','k','LineWidth',2,'Normalization','cdf');
    histogram(dp_VS,-4:histbinsize:4,'DisplayStyle','stairs',...
        'EdgeColor',colors(irate,:),'LineWidth',2,'Normalization','cdf');
    
    xlim([-2 2])
    ylim([0 1])
    xlabel('d prime')
    axis square
    title([num2str(AMrates(irate)) ' Hz'])
end

suptitle('d'' by AM rate, colors: high RS')

print_eps_kp(gcf,fullfile(fn.figs,'MPHclass','dprime_AMrate_highRS'))


end



%%  Compare performance of classifiers

keyboard

dp_L1o = [];
dp_t10 = [];

for iUn = 1:size(Data,1)
    
    for irate = 1:size(Data,2)
        dp_L1o  =  [dp_L1o; Res(iUn,irate).L1o.dprime(:,2)];
        dp_t10  =  [dp_t10; Res(iUn,irate).t10.dprime(:,2)];
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


