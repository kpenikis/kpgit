function quickPlotPopNormSpeech(RERUN)
% 
% plotNormMPHs
%
%  Intended to help categorize response types. Based on plotNormMPHs. 
%
%  Later extract true segments of speech akin to MPH. For now, plot just
%  z-scored PSTH of a few segments of speech responses. 
% 
%  Sort units using the zThLatAbv method.
% 
% 

global AMrates 


Duration  = 500;
Stimuli   = 1:6;

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
% spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

savedir = fullfile(fn.figs,'PopulationSpeechSegments');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
if RERUN
    
    % First X ms of each trial, no folding
    FR_vec    = nan(numel(UnitData),Duration,length(Stimuli));
    zFR_vec   = nan(numel(UnitData),Duration,length(Stimuli));
    
    for iUn = 1:numel(UnitData)
                
        % - - - -   Collect data   - - - - -
        
        get_trial_data_speech
        
        FR_vec(iUn,:,:)  = iu_FRvec;
        zFR_vec(iUn,:,:) = iu_zFRvec;
        
        
    end %iUn
    
    % Save MPH data
    save(fullfile(savedir,'OnsetData'),'zFR_vec','FR_vec','RMS','-v7.3')
    
else
    
    
    
    %% Instead of re-running, load saved MPH data
    
    load(fullfile(savedir,'OnsetData'))
    
    % Load an Info file too
    filename = sprintf( '%s_sess-%s_Info',     UnitData(1).Subject,UnitData(1).Session); load(fullfile(fn.processed,UnitData(1).Subject,filename));
    
end



%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];


%%

SortBy =  'latency'; 'zThLatAbv'; 'signsum'; 'mean'; 'range'; 
thresh = 0.2;
inWin  = 0.25;
zBoundaries = [-1 0 1 2];

%~~~~~~~~~~~~
%   Fig 1
%~~~~~~~~~~~~

hf1 = figure;
set(gcf,'Position',fullscreen)
hold on

for ist=1:size(zFR_vec,3)
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(zFR_vec(:,:,ist),2));
        case 'range'
            [~,i_sorted] = sort(range(zFR_vec(:,:,ist),2));
        case 'signsum'
            [~,i_sorted] = sort(sum(sign(zFR_vec(:,:,ist)),2)) ; 
        case 'zThLatAbv'
            respdata = [];
            for iUn=1:size(zFR_vec,1)
                clear uThZ uLat uAbv
                if all(isnan(zFR_vec(iUn,:,ist)))
                    respdata = [respdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = zBoundaries(find(max(zFR_vec(iUn,:,ist),[],2) > zBoundaries,1,'last'));
                    %find time exceeded it
                    uLat = find(zFR_vec(iUn,:,ist) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(zFR_vec(iUn,:,ist) > uThZ);
                    
                    % Save unit data
                    respdata = [respdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( respdata, [1 -2 3]);
            
        case 'latency'
            zTh = [-0.5 0.5];
            respdata = [];
            for iUn=1:size(zFR_vec,1)
                clear uThZ uLat uAbv
                if all(isnan(zFR_vec(iUn,:,ist)))
                    respdata = [respdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = zTh(find(max(zFR_vec(iUn,:,ist),[],2) > zTh,1,'last'));
                    %find time exceeded it
                    uLat = find(zFR_vec(iUn,:,ist) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(zFR_vec(iUn,:,ist) > uThZ);
                    
                    % Save unit data
                    respdata = [respdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( respdata, [1 -2]);
    end
    
    % Now plot
    subplot(4,size(zFR_vec,3),ist);
    plot(RMS(ist,:),'k','LineWidth',2)
    xlim([0 Duration])
    ylim([0 0.04])
    set(gca,'ytick',[],'xtick',[],'Color','none')
    title(Info.stim_ID_key{ist})
    box off
    
    subplot(4,size(zFR_vec,3),ist+size(zFR_vec,3).*[1 2 3]);
    
    zdata = zFR_vec(i_sorted,:,ist) ;
    
    ndp = sum(sum(isnan(zdata),2)==0);
    
    imagesc(zdata(1:ndp,:))
%     caxis([-1 1])
    caxis([min(zBoundaries) max(zBoundaries)])
%     caxis(max(abs(zBoundaries)).*[-1 1])
    cmocean('balance','pivot',0)
%     colorbar
        
    xlim([0 Duration])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',0:250:Duration,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    if strcmp(SortBy,'zThLatAbv')
        [r,c]=find(cumsum(respdata(i_sorted,1)==zBoundaries)==1);
        set(gca,'ytick',r,'yticklabel',zBoundaries)
    end
    if ist==1
        ylabel('Unit #   |   z-scored FR, norm within each unit')
        xlabel('Time from stim onset')
    end
end %ist 


% colorbar('east')


% Save
print_eps_kp(hf1,fullfile(savedir,['zOnset_' SortBy '' ]))



end





