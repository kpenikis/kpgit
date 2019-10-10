function plotPopNormMPH_v1(RERUN)
% 
% plotNormMPHs
%
%  Original plots from PopulationVariability, now alone and with options
%  for sorting the units.
% 
%  Intended to help categorize response types.
% 

global AMrates 


alfaVS = 0.001;
trMax  = 40;
        
% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'Units_250'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
% spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load RCorr results
q=load(fullfile(fn.figs,'RCorr','exclOnset','PCMat_10trTemp_b.mat'));
PCMat = q.PCMat;
clear q
if size(PCMat,3) ~= numel(UnitData)
    warning('cant sort by RCorr accuracy')
    keyboard
end

savedir = fullfile(fn.figs,'PopulationTuning');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
if RERUN

FR_vec   = nan(numel(UnitData),500,5);
zFR_vec  = nan(numel(UnitData),500,5);
zFR_Warn = nan(numel(UnitData),500);

for iUn = 1:numel(UnitData)
    
    spkshift = [];
    
    % - - - -   Get raw spike data   - - - - - 
    get_trial_data_posthoc
    
    if isempty(spkshift)
        continue
    end
    
    % - - - -   Pull relevant data   - - - - - 
    
    for stid = 2:6 %find(UnitData(iUn).VSdata_spk(2,:)>13.1)
        
        if isempty(MPH(MPH.ThisStimID==stid,:))
            continue
        end
        thisRaster = vertcat(MPH(MPH.ThisStimID==stid,:).raster{:});
        
        % Normalized FR within this unit, referenced to silence
        zFR_vec(iUn,1:size(thisRaster,2),stid-1) = mean(vertcat(MPH(MPH.ThisStimID==stid,:).zFR{:}),1);
        
        % Get smoothed FR
        FR_vec(iUn,1:size(thisRaster,2),stid-1)  = mean(vertcat(MPH(MPH.ThisStimID==stid,:).FRsmooth{:}),1);
        
        
    end %only significant sync
    
    % Save avg Warn response for this unit
    zFR_Warn(iUn,:) = mean(Warn_zFR,1);
    
end %iUn

% Save MPH data
save(fullfile(savedir,'MPHdata_ownSpkShifts'),'zFR_vec','zFR_Warn','FR_vec','-v7.3')



else
    
%% Instead of re-running, load saved MPH data
load(fullfile(savedir,'MPHdata_ownSpkShifts')) 

end



%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


%%

% Things to sort units by:
%  - # bins above some threshold
%  - mean zFR
%  - RCorr accuracy*    PCMat

SortBy =   'signsum'; 'rcorr'; 'width'; 'EdgeRatio'; 'signsum'; 'srMBE'; 'mean'; 'EdgeRatio'; 'peakwidth'; 'range'; 'rcorr'; 
thresh = 0.2;
inWin  = 0.25;

%~~~~~~~~~~~~
%   Fig 1
%~~~~~~~~~~~~

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
    end
    
    % Now plot
    subplot(1,6,ir+1);
    
    zdata = zFR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir);
    
    ndp = sum(sum(isnan(zdata),2)==0);
    
    imagesc(zdata(1:ndp,:))
    caxis([-1 1])
    cmocean('balance','pivot',0)
%     colorbar
        
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    
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
end

% Now plot
subplot(1,6,1);

zdata = zFR_Warn(i_sorted,:);

ndp = sum(sum(isnan(zdata),2)==0);

imagesc(zdata(1:ndp,:))
caxis([-1 1])
cmocean('balance','pivot',0)
%     colorbar

xlim([0 500])
ylim([0.5 ndp+0.5])
set(gca,'ytick',[],'xtick',[0 500],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
title('Warn')
box off
axis fill

ylabel('Unit #   |   z-scored MPH, norm within each unit')
xlabel('Time in MPH')


% Save
print_eps_kp(hf1,fullfile(savedir,['zMPH_' SortBy '_ownShift']))



keyboard
%%
%~~~~~~~~~~~~
%   Fig 2
%~~~~~~~~~~~~

SortBy =   'mean'; 'width'; 'range'; 'rcorr'; 'signsum'; 
thresh = 0.666;
inWin  = 0.25;

hf2 = figure;
set(gcf,'Position',fullscreen)
hold on

% norm FR

for ir=1:5
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(FR_vec(:,1:ceil(1000/AMrates(ir)),ir),2)); %inWin*ceil...
        case 'width'
            [~,i_sorted] = sort(sum(FR_vec(:,:,ir)>thresh,2));
    end
    
    % Now plot
    subplot(2,5,ir+[0 5]);
    
    fdata = FR_vec(i_sorted,1:ceil(1000/AMrates(ir)),ir);
%     fdata = fdata -min(fdata,[],2);
%     fdata = fdata./max(fdata,[],2);
    
    imagesc(fdata)
%     colorbar
    colormap('bone')
%     caxis([-1 3])
%     cmocean('balance','pivot',0)
        
    xlim([0 ceil(1000/AMrates(ir))])
    ylim([0.5 size(FR_vec,1)+0.5])
    set(gca,'ytick',[],'xtick',[0 ceil(1000/AMrates(ir))],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    if ir==1
        ylabel('Unit #   |   norm FR MPH, norm within each MPH')
        xlabel('Time in MPH')
    end
    title([num2str(AMrates(ir)) ' Hz'])
    box off
    axis fill
    
end


% Save figure
print_eps_kp(hf2,fullfile(savedir,['nMPH_' SortBy]))



end



