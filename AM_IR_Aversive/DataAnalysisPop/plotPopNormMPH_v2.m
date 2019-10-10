function plotPopNormMPH_v2(RERUN)
% 
% plotNormMPHs_v2
%
%  Original plots from PopulationVariability, now separate, with options
%  for sorting the units.
%
%  Now broken off to focus on dynamics of population activity during each
%  AM rate MPH.
% 

global AMrates 

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
widescreen  = [1 scrsz(4) scrsz(3) scrsz(4)/2];

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

alfaVS = 0.001;
trMax  = 40;

zBoundaries = [-1 0 0.25 0.5 1 2];

%%
% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load RCorr results
q=load(fullfile(fn.figs,'RCorr','exclOnset','PCMat_10trTemp_b.mat'));
PCMat = q.PCMat;
clear q
if size(PCMat,3) ~= numel(UnitData)
    warning('cant sort by RCorr accuracy')
    keyboard
end

% savedir = fullfile(fn.figs,'PopulationTuning');
savedir = fullfile(fn.figs,'PopMPH');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
if RERUN
    
    %currently: common spike shift
    
    FR_vec   = nan(numel(UnitData),500,5);
    zFR_vec  = nan(numel(UnitData),500,5);
    FR_Warn  = nan(numel(UnitData),500);
    zFR_Warn = nan(numel(UnitData),500);
    
    sp_trs       = nan(numel(UnitData),500,5,10);
    sp_Warn_trs  = nan(numel(UnitData),500,10);
    
    for iUn = 1:numel(UnitData)
        
        %     spkshift = [];
        
        % - - - -   Get raw spike data   - - - - -
        get_trial_data_posthoc
        
        %     if isempty(spkshift)
        %         continue
        %     end
        
        % - - - -   Pull relevant data   - - - - -
        
        for stid = 2:6 
            
            if isempty(MPH(MPH.ThisStimID==stid,:))
                continue
            end
            
            % Save 10 random trials
            thisRaster = vertcat(MPH(MPH.ThisStimID==stid,:).raster{:});
            ridx = randperm(size(thisRaster,1),10);
            sp_trs(iUn,1:size(thisRaster,2),stid-1,:) = thisRaster(ridx,:)';
            
            % Normalized FR within this unit, referenced to silence
            zFR_vec(iUn,1:size(thisRaster,2),stid-1) = mean(vertcat(MPH(MPH.ThisStimID==stid,:).zFR{:}),1);
            
            % Get smoothed FR
            FR_vec(iUn,1:size(thisRaster,2),stid-1)  = mean(vertcat(MPH(MPH.ThisStimID==stid,:).FRsmooth{:}),1);
            
            
        end %only significant sync
        
        % Save avg Warn response for this unit
        FR_Warn(iUn,:)  = mean(Warn_FR,1);
        zFR_Warn(iUn,:) = mean(Warn_zFR,1);
        
    end %iUn
    
    % Save MPH data
    save(fullfile(savedir,'MPHdata'),'zFR_vec','zFR_Warn','FR_vec','FR_Warn','sp_trs','-v7.3')
    
    
else
    
    % Instead of re-running, load saved MPH data
    
    load(fullfile(savedir,'MPHdata'))
    
end


% - cumulative fr, overlayed
% - psth overlayed

%   repeat with raw spikes, not smoothed *

%   cdf every 10ms 

% - within cell - overlay MPHs 
%   quantify if timing shifts or stays same

% - RATIO of zscore


keyboard

%% 

% figure;
% set(gcf,'Position',widescreen)
% 
% for ir=1:5
%     
%     subplot(1,3,1:2); hold on
%     plot(mean(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',3)
%     plot(mean(zFR_vec(:,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',6)
%     
%     subplot(1,3,3); hold on
%     plot(cumsum(mean(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),1,'omitnan')),'Color',colors(ir+1,:),'LineWidth',3)
%     plot(cumsum(mean(zFR_vec(:,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan')),'Color',colors(ir+1,:),'LineWidth',6)
%     
% end
% 
% subplot(1,3,1:2); hold on
% plot(mean(zFR_Warn,1,'omitnan'),'Color',colors(8,:),'LineWidth',3)
% % plot(mean(zFR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan'),'Color',colors(8,:),'LineWidth',6)
% 
% title('average zFR')
% xlabel('time (ms)')
% 
% subplot(1,3,3); hold on
% plot(cumsum(mean(zFR_Warn,1,'omitnan')),'Color',colors(8,:),'LineWidth',3)
% % plot(cumsum(mean(zFR_Warn(:,1:size(zFR_Warn,2)/2),1,'omitnan')),'Color',colors(8,:),'LineWidth',6)
% 
% title('cumulative zFR')
% xlabel('time (ms)')
% 
% print_eps_kp(gcf,fullfile(savedir,'Pop_zFR_MPH'))
% 
% 
% 
% figure;
% set(gcf,'Position',widescreen)
% 
% for ir=1:5
%     
%     subplot(1,3,1:2); hold on
%     plot(mean(FR_vec(:,1:ceil(1000/AMrates(ir)),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',3)
%     plot(mean(FR_vec(:,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',6)
%     
%     subplot(1,3,3); hold on
%     plot(cumsum(sum(FR_vec(:,1:ceil(1000/AMrates(ir)),ir),1,'omitnan')./1000),'Color',colors(ir+1,:),'LineWidth',3)
%     plot(cumsum(sum(FR_vec(:,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan')./1000),'Color',colors(ir+1,:),'LineWidth',6)
%     
% end
% 
% subplot(1,3,1:2); hold on
% plot(mean(FR_Warn,1,'omitnan'),'Color',colors(8,:),'LineWidth',3)
% % plot(mean(FR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan'),'Color',colors(8,:),'LineWidth',6)
% 
% title('average FR')
% xlabel('time (ms)')
% 
% subplot(1,3,3); hold on
% plot(cumsum(sum(FR_Warn,1,'omitnan')./1000),'Color',colors(8,:),'LineWidth',3)
% % plot(cumsum(sum(FR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan')./1000),'Color',colors(8,:),'LineWidth',6)
% 
% title('cumulative # spikes')
% xlabel('time (ms)')
% 
% print_eps_kp(gcf,fullfile(savedir,'Pop_FR_MPH'))


%% 

% Individual cells  -  latency shifts 
% 
% 

for iUn = 1:size(zFR_vec,1)
    
    if max(zFR_vec(iUn,1:ceil(1000/AMrates(2)),2))<1
        continue
    end
    
    hf1=figure;
    set(gcf,'Position',widescreen)

for ir=1:5
    
    subplot(1,3,1:2); hold on
    plot(mean(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',3)
    plot(mean(zFR_vec(iUn,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',6)
    
    subplot(1,3,3); hold on
    plot(cumsum(mean(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),1,'omitnan')),'Color',colors(ir+1,:),'LineWidth',3)
    plot(cumsum(mean(zFR_vec(iUn,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan')),'Color',colors(ir+1,:),'LineWidth',6)
    
end

subplot(1,3,1:2); hold on
plot(mean(zFR_Warn(iUn,:),1,'omitnan'),'Color',colors(8,:),'LineWidth',3)
% plot(mean(zFR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan'),'Color',colors(8,:),'LineWidth',6)

title('average zFR')
xlabel('time (ms)')

subplot(1,3,3); hold on
plot(cumsum(mean(zFR_Warn(iUn,:),1,'omitnan')),'Color',colors(8,:),'LineWidth',3)
% plot(cumsum(mean(zFR_Warn(:,1:size(zFR_Warn,2)/2),1,'omitnan')),'Color',colors(8,:),'LineWidth',6)

title('cumulative zFR')
xlabel('time (ms)')




hf2=figure;
set(gcf,'Position',widescreen)

for ir=1:5
    
    subplot(1,3,1:2); hold on
    plot(mean(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',3)
    plot(mean(FR_vec(iUn,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan'),'Color',colors(ir+1,:),'LineWidth',6)
    
    subplot(1,3,3); hold on
    plot(cumsum(sum(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),1,'omitnan')./1000),'Color',colors(ir+1,:),'LineWidth',3)
    plot(cumsum(sum(FR_vec(iUn,1:ceil(1000/AMrates(ir)/2),ir),1,'omitnan')./1000),'Color',colors(ir+1,:),'LineWidth',6)
    
end

subplot(1,3,1:2); hold on
plot(mean(FR_Warn(iUn,:),1,'omitnan'),'Color',colors(8,:),'LineWidth',3)
% plot(mean(FR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan'),'Color',colors(8,:),'LineWidth',6)

title('average FR')
xlabel('time (ms)')

subplot(1,3,3); hold on
plot(cumsum(sum(FR_Warn(iUn,:),1,'omitnan')./1000),'Color',colors(8,:),'LineWidth',3)
% plot(cumsum(sum(FR_Warn(:,1:size(FR_Warn,2)/2),1,'omitnan')./1000),'Color',colors(8,:),'LineWidth',6)

title('cumulative # spikes')
xlabel('time (ms)')


% print_eps_kp(hf1,fullfile(savedir,'IndUnits',['Pop_zFR_MPH_' num2str(iUn)]))
% print_eps_kp(hf2,fullfile(savedir,'IndUnits',['Pop_FR_MPH_'  num2str(iUn)]))

keyboard

end

%% Rank correlation of sorting across AM rates

Sorted = nan(size(zFR_vec,1),5);

for ir=1:5
    
    % Sort here
    respdata = [];
    
    % do it the slow but easy way
    for iUn=1:size(zFR_vec,1)
        
        clear uThZ uLat uAbv
        
        if all(isnan(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir)))
            respdata = [respdata; nan nan nan ];
        else
            
            %find highest thresh exceeded
            uThZ = zBoundaries(find(max(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),[],2) > zBoundaries,1,'last'));
            %find time exceeded it
            uLat = find(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uThZ,1,'first');
            %find ms above
            uAbv = sum(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uThZ);
            
            % Save unit data
            respdata = [respdata; uThZ uLat uAbv ];
        end
    end
    [~,i_sorted] = sortrows( respdata, [1 -2 -3]);
    [i_sorted_new,Sorted(:,ir)] = sort(flipud(i_sorted));
    
    Sorted(respdata(:,1)<0.5,ir) = nan;
    
end

[r,p] = corr(Sorted,'type','Spearman')


[ss,iss] = sort(i_sorted)
flipud(Sorted)



end





%%  -----    zFR    -----
function WhichPlot1

SortBy =   'zThLatAbv'; 'signsum'; 'rcorr'; 'width'; 'EdgeRatio'; 'signsum'; 'srMBE'; 'mean'; 'EdgeRatio'; 'peakwidth'; 'range'; 'rcorr'; 

zBoundaries = [-1 0 0.25 0.5 1 2];

thresh = 0.2;
inWin  = 0.25;

hf1 = figure;
set(gcf,'Position',fullscreen)
hold on

for ir=1:5
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),2));
        case 'range'
            [~,i_sorted] = sort(range(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir),2));
        case 'width'
            [~,i_sorted] = sort(sum(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir)>thresh,2));
        case 'signsum'
            [~,i_sorted] = sort(sum(sign(zFR_vec(:,1:ceil(1000/AMrates(ir)),ir)),2)) ; 
        case 'peakwidth'
            pw = zeros(size(zFR_vec,1),1);
            for iUn=1:size(zFR_vec,1)
                ii = find(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir)>thresh);
                if ~isempty(ii)
                    pw(iUn) = sum(diff(ii));
                end
            end
            [~,i_sorted] = sort(pw);
        case 'EdgeRatio'
            edg = mean(zFR_vec(:,[1:ceil(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):end],ir),2);
            mid = mean(zFR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),2);
            [~,i_sorted] = sort((mid-edg)./(mid+edg));
        case 'srMBE'
            mid = round(mean(zFR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),2),1);
            beg = round(mean(zFR_vec(:,1:ceil(0.25*1000/AMrates(ir)),ir),2),1);
            nnd = round(mean(zFR_vec(:,ceil(0.75*1000/AMrates(ir)):end,ir),2),1);
%             [~,i_sorted] = sortrows([mid-mean([beg nnd],2) mid beg nnd],[1 2 3 4]);
            [~,i_sorted] = sortrows([mid mean([beg nnd],2)],[1 2]);
        case 'rcorr'
            [~,i_sorted] = sort(permute(PCMat(ir+1,ir+1,:),[3 1 2]));
            
        case 'zThLatAbv'
            respdata = [];
            uAbvRnd = ceil(1000/AMrates(ir)/2); %1/4 of period
            
            % do it the slow but easy way
            for iUn=1:size(zFR_vec,1)
                
                clear uThZ uLat uAbv
                
                if all(isnan(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir)))
                    respdata = [respdata; nan nan nan ];
                else
                    
                    %find highest thresh exceeded
                    uThZ = zBoundaries(find(max(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),[],2) > zBoundaries,1,'last'));
                    %find time exceeded it
                    uLat = find(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uThZ);
%                     uAbv = ceil(sum(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > (uThZ-0.5))/uAbvRnd)*uAbvRnd;
                    
                    % Save unit data
                    respdata = [respdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( respdata, [1 -2 -3]);
    end 
    
    % Sort data
    zdata = [zFR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir) zFR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir)];
    
    ndp = sum(sum(isnan(zdata),2)==0);
    
    % Label NS
    flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
    
    % Now plot
    subplot(1,6,ir+1);
    
    imagesc(zdata(1:ndp,:))
%     caxis([-1 1])
    caxis([min(zBoundaries) max(zBoundaries)])
%     caxis(max(abs(zBoundaries)).*[-1 1])
    cmocean('balance','pivot',0)
    
    % Add markers to label NS cells
    hold on
    plot(0,find(flagNS),'.g')
    plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.g')
    
    xlim([0 2*ceil(1000/AMrates(ir))])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    if strcmp(SortBy,'zThLatAbv')
        [r,c]=find(cumsum(respdata(i_sorted,1)==zBoundaries)==1);
        set(gca,'ytick',r,'yticklabel',zBoundaries)
    end
    
end


%~~~~  Warn response  ~~~~%

% Sort here
switch SortBy
    case 'mean'
        [~,i_sorted] = sort(mean(zFR_Warn,2));
    case 'range'
        [~,i_sorted] = sort(range(zFR_Warn,2));
    case 'width'
        [~,i_sorted] = sort(sum(zFR_Warn>thresh,2));
    case 'signsum'
        [~,i_sorted] = sort(sum(sign(zFR_Warn),2)) ;
    case 'peakwidth'
        pw = zeros(size(zFR_Warn,1),1);
        for iUn=1:size(zFR_Warn,1)
            ii = find(zFR_Warn(iUn,1:ceil(1000/AMrates(ir)),ir)>thresh);
            if ~isempty(ii)
                pw(iUn) = sum(diff(ii));
            end
        end
        [~,i_sorted] = sort(pw);
    case 'EdgeRatio'
        edg = mean(zFR_Warn(:,[ 1:125 376:end ]),2);
        mid = mean(zFR_Warn(:,126:375),2);
        [~,i_sorted] = sort((mid-edg)./(mid+edg));
    case 'srMBE'
        mid = round(mean(zFR_Warn(:,126:375),2),1);
        beg = round(mean(zFR_Warn(:,  1:125),2),1);
        nnd = round(mean(zFR_Warn(:,376:end),2),1);
        %             [~,i_sorted] = sortrows([mid-mean([beg nnd],2) mid beg nnd],[1 2 3 4]);
        [~,i_sorted] = sortrows([mid mean([beg nnd],2)],[1 2]);
    case 'rcorr'
        [~,i_sorted] = sort(permute(PCMat(1,1,:),[3 1 2]));
        
    case 'zThLatAbv'
        respdata = [];
        uAbvRnd = ceil(500/2); %1/4 of period
        % do it the slow but easy way
        for iUn=1:size(zFR_Warn,1)
            clear uThZ uLat uAbv
            if all(isnan(zFR_Warn(iUn,:)))
                respdata = [respdata; nan nan nan ];
            else
                %find highest thresh exceeded
                uThZ = zBoundaries(find(max(zFR_Warn(iUn,:),[],2) > zBoundaries,1,'last'));
                %find time exceeded it
                uLat = find(zFR_Warn(iUn,:) > uThZ,1,'first');
                %find ms above
                uAbv = sum(zFR_Warn(iUn,:) > uThZ);
%                 uAbv = ceil(sum(zFR_Warn(iUn,:) > (uThZ-0.5))/uAbvRnd)*uAbvRnd;
                
                % Save unit data
                respdata = [respdata; uThZ uLat uAbv ];
            end
        end
        [~,i_sorted] = sortrows( respdata, [1 -2 3]);
end

% Sort data
zdata = zFR_Warn(i_sorted,:);

ndp = sum(sum(isnan(zdata),2)==0);

% Label NS
flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;

% Now plot
subplot(1,6,1);

imagesc(zdata(1:ndp,:))
% caxis(max(abs(zBoundaries)).*[-1 1])
caxis([min(zBoundaries) max(zBoundaries)])
cmocean('balance','pivot',0)

% Add markers to label NS cells
hold on
plot(0,find(flagNS),'.g')
plot(500,find(flagNS),'.g')

xlim([0 500])
ylim([0.5 ndp+0.5])
set(gca,'ytick',[],'xtick',[0 500],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
title('Warn')
box off
axis fill
if strcmp(SortBy,'zThLatAbv')
    [r,c]=find(cumsum(respdata(i_sorted,1)==zBoundaries)==1);
    set(gca,'ytick',r,'yticklabel',zBoundaries)
end

ylabel('Unit #   |   z-scored MPH, norm within each unit')
xlabel('Time in MPH')

% colorbar('east')

% Save
print_eps_kp(hf1,fullfile(savedir,['PopMPH_z_' SortBy '' ]))

end

%%  -----   log FR   -----
function WhichPlot2

SortBy =   'ThLatAbv'; 'signsum'; 'rcorr'; 'width'; 'EdgeRatio'; 'signsum'; 'srMBE'; 'mean'; 'EdgeRatio'; 'peakwidth'; 'range'; 'rcorr'; 

frBoundaries = [-1 round(10.^[0.5 1 1.25 1.5 1.75]) ];

thresh = 0.2;
inWin  = 0.25;

hf2 = figure;
set(gcf,'Position',fullscreen)
hold on

for ir=1:5
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(FR_vec(:,1:ceil(1000/AMrates(ir)),ir),2));
        case 'range'
            [~,i_sorted] = sort(range(FR_vec(:,1:ceil(1000/AMrates(ir)),ir),2));
        case 'width'
            [~,i_sorted] = sort(sum(FR_vec(:,1:ceil(1000/AMrates(ir)),ir)>thresh,2));
        case 'signsum'
            [~,i_sorted] = sort(sum(sign(FR_vec(:,1:ceil(1000/AMrates(ir)),ir)),2)) ; 
        case 'peakwidth'
            pw = zeros(size(FR_vec,1),1);
            for iUn=1:size(FR_vec,1)
                ii = find(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir)>thresh);
                if ~isempty(ii)
                    pw(iUn) = sum(diff(ii));
                end
            end
            [~,i_sorted] = sort(pw);
        case 'EdgeRatio'
            edg = mean(FR_vec(:,[1:ceil(0.25*1000/AMrates(ir)) ceil(0.75*1000/AMrates(ir)):end],ir),2);
            mid = mean(FR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),2);
            [~,i_sorted] = sort((mid-edg)./(mid+edg));
        case 'srMBE'
            mid = round(mean(FR_vec(:,ceil(0.25*1000/AMrates(ir)):ceil(0.75*1000/AMrates(ir)),ir),2),1);
            beg = round(mean(FR_vec(:,1:ceil(0.25*1000/AMrates(ir)),ir),2),1);
            nnd = round(mean(FR_vec(:,ceil(0.75*1000/AMrates(ir)):end,ir),2),1);
%             [~,i_sorted] = sortrows([mid-mean([beg nnd],2) mid beg nnd],[1 2 3 4]);
            [~,i_sorted] = sortrows([mid mean([beg nnd],2)],[1 2]);
        case 'rcorr'
            [~,i_sorted] = sort(permute(PCMat(ir+1,ir+1,:),[3 1 2]));
            
        case 'ThLatAbv'
            respdata = [];
            uAbvRnd = ceil(1000/AMrates(ir)/2); %1/4 of period
            
            % do it the slow but easy way
            for iUn=1:size(FR_vec,1)
                
                clear uTh uLat uAbv
                
                if all(isnan(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir)))
                    respdata = [respdata; nan nan nan ];
                else
                    
                    %find highest thresh exceeded
                    uTh = frBoundaries(find(max(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir),[],2) > frBoundaries,1,'last'));
                    %find time exceeded it
                    uLat = find(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uTh,1,'first');
                    %find ms above
                    uAbv = sum(FR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > uTh);
%                     uAbv = ceil(sum(zFR_vec(iUn,1:ceil(1000/AMrates(ir)),ir) > (uThZ-0.5))/uAbvRnd)*uAbvRnd;
                    
                    % Save unit data
                    respdata = [respdata; uTh uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( respdata, [1 -2 -3]);
    end 
    
    
    
    % Sort data
    frdata = [FR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir) FR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir)];
    
    ndp = sum(sum(isnan(frdata),2)==0);
    
    % Label NS
    flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
    
    % Now plot
    subplot(1,6,ir+1);
    
    imagesc(log10(frdata(1:ndp,:)))
%     caxis([-1 1])
    caxis([0 log10(max(frBoundaries))])
%     caxis(max(abs(zBoundaries)).*[-1 1])
    cmocean('gray')
    
    % Add markers to label NS cells
    hold on
    plot(0,find(flagNS),'.g')
    plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.g')
    
    xlim([0 2*ceil(1000/AMrates(ir))])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir)) 2*ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    if strcmp(SortBy,'ThLatAbv')
        [r,c]=find(cumsum(respdata(i_sorted,1)==frBoundaries)==1);
        set(gca,'ytick',r,'yticklabel',frBoundaries)
    end
    
end


%~~~~  Warn response  ~~~~%

% Sort here
switch SortBy
    case 'mean'
        [~,i_sorted] = sort(mean(FR_Warn,2));
    case 'range'
        [~,i_sorted] = sort(range(FR_Warn,2));
    case 'width'
        [~,i_sorted] = sort(sum(FR_Warn>thresh,2));
    case 'signsum'
        [~,i_sorted] = sort(sum(sign(FR_Warn),2)) ;
    case 'peakwidth'
        pw = zeros(size(FR_Warn,1),1);
        for iUn=1:size(FR_Warn,1)
            ii = find(FR_Warn(iUn,1:ceil(1000/AMrates(ir)),ir)>thresh);
            if ~isempty(ii)
                pw(iUn) = sum(diff(ii));
            end
        end
        [~,i_sorted] = sort(pw);
    case 'EdgeRatio'
        edg = mean(FR_Warn(:,[ 1:125 376:end ]),2);
        mid = mean(FR_Warn(:,126:375),2);
        [~,i_sorted] = sort((mid-edg)./(mid+edg));
    case 'srMBE'
        mid = round(mean(FR_Warn(:,126:375),2),1);
        beg = round(mean(FR_Warn(:,  1:125),2),1);
        nnd = round(mean(FR_Warn(:,376:end),2),1);
        %             [~,i_sorted] = sortrows([mid-mean([beg nnd],2) mid beg nnd],[1 2 3 4]);
        [~,i_sorted] = sortrows([mid mean([beg nnd],2)],[1 2]);
    case 'rcorr'
        [~,i_sorted] = sort(permute(PCMat(1,1,:),[3 1 2]));
        
    case 'ThLatAbv'
        respdata = [];
        uAbvRnd = ceil(500/2); %1/4 of period
        % do it the slow but easy way
        for iUn=1:size(FR_Warn,1)
            clear uTh uLat uAbv
            if all(isnan(FR_Warn(iUn,:)))
                respdata = [respdata; nan nan nan ];
            else
                %find highest thresh exceeded
                uTh = frBoundaries(find(max(FR_Warn(iUn,:),[],2) > frBoundaries,1,'last'));
                %find time exceeded it
                uLat = find(FR_Warn(iUn,:) > uTh,1,'first');
                %find ms above
                uAbv = sum(FR_Warn(iUn,:) > uTh);
%                 uAbv = ceil(sum(zFR_Warn(iUn,:) > (uThZ-0.5))/uAbvRnd)*uAbvRnd;
                
                % Save unit data
                respdata = [respdata; uTh uLat uAbv ];
            end
        end
        [~,i_sorted] = sortrows( respdata, [1 -2 -3]);
end



% Sort data
frdata = FR_Warn(i_sorted,:);

ndp = sum(sum(isnan(frdata),2)==0);

% Label NS
flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;

% Now plot
subplot(1,6,1);

imagesc(log10(frdata(1:ndp,:)))
caxis([0 log10(max(frBoundaries))])
cmocean('gray')

% Add markers to label NS cells
hold on
plot(0,find(flagNS),'.g')
plot(500,find(flagNS),'.g')

xlim([0 500])
ylim([0.5 ndp+0.5])
set(gca,'ytick',[],'xtick',[0 500],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
title('Warn')
box off
axis fill
if strcmp(SortBy,'ThLatAbv')
    [r,c]=find(cumsum(respdata(i_sorted,1)==frBoundaries)==1);
    set(gca,'ytick',r,'yticklabel',frBoundaries)
end

ylabel('Unit #   |   MPH, log FR')
xlabel('Time in MPH')

% colorbar('east')

% Save
print_eps_kp(hf2,fullfile(savedir,['PopMPH_log_' SortBy '' ]))


end
