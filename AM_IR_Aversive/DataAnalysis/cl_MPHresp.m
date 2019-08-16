function cl_MPHresp
% 


global AMrates trMin rateVec_AC rateVec_DB

alfaVS = 0.001;
trMin  = 5;
trMax  = 40;
        
% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load classifier results
load(fullfile(fn.processed,'MPHclassifier','ClassData'));

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;
clear q

AMrates = [2 4 8 16 32];


%%

zFR_vec  = zeros(numel(UnitData),500,5);
pkStats  = nan(numel(UnitData),3,5);

for iUn = 1:numel(UnitData)
    
    % Get raw spike data   
    get_trial_data_posthoc
    
    
    % Get zscored MPH
    for stid = find(UnitData(iUn).VSdata_spk(2,:)>13.1)
        
        thisRaster = vertcat(MPH(MPH.ThisStimID==stid,:).raster{:});
        
        % Only analyze trials with spikes
        trials = find(sum(thisRaster,2)>0)';
        
        if isempty(trials)
            continue
        end
        
        thisMPH = mean(vertcat(MPH(MPH.ThisStimID==stid,:).zFR{:}),1);
        zFR_vec(iUn,1:size(thisRaster,2),stid-1) = thisMPH;
        
        if max(zFR_vec(iUn,1:size(thisRaster,2),stid-1)) < 0.5
            continue
        end
        
        pkheight=[]; pklatency=[]; pkwidth=[];
        [pkheight,pklatency,pkwidth]= findpeaks(zFR_vec(iUn,1:size(thisRaster,2),stid-1),'MinPeakHeight',0.5,'MinPeakProminence',0.5,'WidthReference','halfprom');
        
        if ~isempty(pkheight)
            
            threshold = (max(thisMPH)-min(thisMPH))/2 + min(thisMPH);
            excPeak = find(thisMPH>threshold);
            
            pkStats(iUn,:,stid-1) = [max(pkheight) mean(excPeak) max(excPeak)-min(excPeak)];
            
%             figure(4); clf;
%             plot(thisMPH,'k')
%             hold on
%             plot([min(excPeak) max(excPeak)],[threshold threshold],'g')
%             plot(mean(excPeak),max(pkheight),'*g')
%             pause(2)
        end
        
    end %only significant sync
    
end %iUn

save(fullfile(fn.processed,'MPHclassifier','excPeakData'),'pkStats','-v7.3')

%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
twothirds   = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)];
widescreen  = [1 scrsz(4) scrsz(3) scrsz(4)/3*2];
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


%%
% pkStats = [  z-height   latency   width  ]

%~~~~~~~~~~~~
% various plots to look for relationship of d' to MPH response type 


% latency vs width of MPH peak, colored by d', only dps with a moderately large response in this MPH
hf1 = figure; 
set(gcf,'Position',fullscreen)

for ir = 1:5
    
    subplot(2,3,ir)
    hold on
    title([num2str(AMrates(ir)) ' Hz'])
    
    for iUn = 1:numel(UnitData)
        
%         if pkStats(iUn,1,ir) > 1
            scatter( pkStats(iUn,2,ir), pkStats(iUn,3,ir), 80, max(Data(iUn,ir).Res_L1o.dprime(:,2)), 'LineWidth',2 )
%         end
        
    end
    xlabel('MPH peak latency')
    ylabel('MPH peak width')
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0 ceil(1000/AMrates(ir))])
    cmocean('phase')
    set(gca,'clim',[0 2])
    axis square
end
subplot(2,3,6)
colorbar
cmocean('phase')
set(gca,'clim',[0 2])


% latency vs width of MPH peak, colored by d', only dps > X
% latency or width vs d', only dps > x
hf2 = figure;
set(gcf,'Position',fullscreen)

for ir = 1:5
    
    subplot(2,3,ir)
    hold on
    title([num2str(AMrates(ir)) ' Hz'])
    
    for iUn = 1:numel(UnitData)
        
%         if ~isempty(Data(iUn,ir).Res_L1o) && any(Data(iUn,ir).Res_L1o.dprime(:,2)>0.5)
%             scatter(pkStats(iUn,2,ir),pkStats(iUn,3,ir),80,max(Data(iUn,ir).Res_L1o.dprime(:,2)),'LineWidth',2)
%         end
        if pkStats(iUn,1,ir) > 1
            plot(pkStats(iUn,3,ir), max(Data(iUn,ir).Res_L1o.dprime(:,2)), '*k')
        end
        
    end
    xlabel('MPH peak width')
    xlim([0 ceil(1000/AMrates(ir))])
    ylabel('best d''')
    ylim([0 2])
%     ylabel('width')
%     ylim([0 ceil(1000/AMrates(ir))])
%     cmocean('phase')
%     set(gca,'clim',[0 6])
    axis square
end

subplot(2,3,6)
colorbar
cmocean('phase')
set(gca,'clim',[0 6])



%~~~~~~~~~~~~
% MPH plots sorted by max dprime

dp_edges = 0:0.75:4.5;

for ir = 1:5
    
    hf(ir) = figure;
    set(gcf,'Position',fullscreen)
    
    for iUn = 1:numel(UnitData)
        
        if ~isempty(Data(iUn,ir).Res_L1o) && ~all(isnan(Data(iUn,ir).Res_L1o.dprime(:,2)))
            
            dpsp = find((dp_edges-max(Data(iUn,ir).Res_L1o.dprime(:,2)))>0,1,'first');
                        
            subplot(2,4,dpsp)
            hold on
            plot(zFR_vec(iUn,:,ir),'k')
            title(sprintf('dp < %0.2f',dp_edges(dpsp)))
            axis square
            xlim([0 ceil(1000/AMrates(ir))])
            ylim([-1 10])
            
        end
    end
    suptitle([num2str(AMrates(ir)) ' Hz z-scored MPHs'])
end


%~~~~~~~~~~~~
% visualize tuning curves with dprimes

for iUn = 1:20:numel(UnitData)
    
    dps = nan(1,5);
    
    for ir = 1:5
        
        if ~isempty(Data(iUn,ir).Res_L1o) && ~all(isnan(Data(iUn,ir).Res_L1o.dprime(:,2)))
            
            dps(ir) = max(Data(iUn,ir).Res_L1o.dprime(:,2));
            
        end
    end
    
    figure;
    plot(1:5,dps,'b','LineWidth',2)
    hold on
    plot(1:5,UnitData(iUn).VSdata_spk(2,2:6) ./max(UnitData(iUn).VSdata_spk(2,2:6)),'r')
    plot(1:5,UnitData(iUn).FR_nrm(2:6),'k')
    ylim([-0.5 max([dps 1])])
    title(num2str(iUn))
end




%~~~~~~~~~~~~
% no correlation between BMFs, or with AM rate of max dprime

BMF_FR = nan(numel(UnitData),1);
BMF_VS = nan(numel(UnitData),1);
max_dp = nan(numel(UnitData),2);

for iUn = 1:numel(UnitData)
    
    dps = nan(1,5);
    
    for ir = 1:5
        
        if ~isempty(Data(iUn,ir).Res_L1o) && ~all(isnan(Data(iUn,ir).Res_L1o.dprime(:,2)))
            
            dps(ir) = median(Data(iUn,ir).Res_L1o.dprime(:,2)); % nor max, nor mean
            
        end
        
        [mdp,idp] = max(dps);
        max_dp(iUn,:) = [mdp idp];
        
        if ~isempty(UnitData(iUn).iBMF_FR)
            BMF_FR(iUn) = UnitData(iUn).iBMF_FR;
        end
        if ~isempty(UnitData(iUn).iBMF_VS)
            BMF_VS(iUn) = UnitData(iUn).iBMF_VS;
        end
        
    end %ir
end %iUn

figure; 
plot(BMF_FR+0.4*rand(numel(UnitData),1)-0.2,max_dp(:,2)+0.4*rand(numel(UnitData),1)-0.2,'.k')
xlim([0 6])
ylim([0 6])
xlabel('BMF FR')
ylabel('rate of max d''')

figure; 
plot(BMF_VS+0.4*rand(numel(UnitData),1)-0.2,max_dp(:,2)+0.4*rand(numel(UnitData),1)-0.2,'.r')
xlim([0 6])
ylim([0 6])
xlabel('BMF VS')
ylabel('rate of max d''')

figure; 
plot(BMF_VS+0.4*rand(numel(UnitData),1)-0.2,BMF_FR+0.4*rand(numel(UnitData),1)-0.2,'.g')
xlim([0 6])
ylim([0 6])
xlabel('BMF VS')
ylabel('BMF FR')


keyboard


% Save figures

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end
print_eps_kp(hf1,fullfile(savedir,'VariabilitySummary_rate'))
print_eps_kp(hf2,fullfile(savedir,'VariabilitySummary_meanphase'))
print_eps_kp(hf3,fullfile(savedir,'VariabilitySummary_zscMPH'))


end





