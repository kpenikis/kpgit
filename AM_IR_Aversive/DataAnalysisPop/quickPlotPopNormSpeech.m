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
useFR = 'log';
switch useFR
    case 'z'
        ThisData    = zFR_vec;
        Boundaries  = [-1 0 0.25 0.5 1 2];
    case 'log'
        ThisData    = FR_vec;
        Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
end

SortBy =  'zThLatAbv'; 'latency'; 'signsum'; 'mean'; 'range'; 
thresh = 0.2;
inWin  = 0.25;

%~~~~~~~~~~~~
%   Fig 1
%~~~~~~~~~~~~

hf1 = figure;
set(gcf,'Position',fullscreen)
hold on

for ist=1:size(ThisData,3)
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(ThisData(:,:,ist),2));
        case 'range'
            [~,i_sorted] = sort(range(ThisData(:,:,ist),2));
        case 'signsum'
            [~,i_sorted] = sort(sum(sign(ThisData(:,:,ist)),2)) ; 
        case 'zThLatAbv'
            sortdata = [];
            for iUn=1:size(ThisData,1)
                clear uThZ uLat uAbv
                if all(isnan(ThisData(iUn,:,ist)))
                    sortdata = [sortdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = Boundaries(find(max(ThisData(iUn,:,ist),[],2) > Boundaries,1,'last'));
                    %find time exceeded it
                    uLat = find(ThisData(iUn,:,ist) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(ThisData(iUn,:,ist) > uThZ);
                    
                    % Save unit data
                    sortdata = [sortdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( sortdata, [1 -2 3]);
            
        case 'latency'
            zTh = [-0.5 0.5];
            sortdata = [];
            for iUn=1:size(ThisData,1)
                clear uThZ uLat uAbv
                if all(isnan(ThisData(iUn,:,ist)))
                    sortdata = [sortdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = zTh(find(max(ThisData(iUn,:,ist),[],2) > zTh,1,'last'));
                    %find time exceeded it
                    uLat = find(ThisData(iUn,:,ist) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(ThisData(iUn,:,ist) > uThZ);
                    
                    % Save unit data
                    sortdata = [sortdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( sortdata, [1 -2]);
    end
    
    
    % Plot stim RMS
    subplot(4,size(ThisData,3),ist);
    plot(RMS(ist,:),'k','LineWidth',2)
    xlim([0 Duration])
    ylim([0 0.04])
    set(gca,'ytick',[],'xtick',[],'Color','none')
    title(Info.stim_ID_key{ist})
    box off
    
    % Plot responses
    subplot(4,size(ThisData,3),ist+size(ThisData,3).*[1 2 3]);
    
    plotdata = ThisData(i_sorted,:,ist) ;
    
    ndp = sum(sum(isnan(plotdata),2)==0);
    
    % Render plot
    switch useFR
        case 'z'
            imagesc(plotdata(1:ndp,:))
            caxis([min(Boundaries) max(Boundaries)])
            cmocean('balance','pivot',0) %curl
        case 'log'
            imagesc(log10(plotdata(1:ndp,:)))
            caxis([0 log10(max(Boundaries))+0.25])
            cmocean('gray')
    end
    
    xlim([0 Duration])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',0:250:Duration,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
    box off
    axis fill
    if strcmp(SortBy,'zThLatAbv')
        BoundMat = cumsum(sortdata(i_sorted(1:ndp),1)>=Boundaries);
        ytplc = []; ytlab = [];
        for ib = numel(Boundaries):-1:1
            yUn = find(BoundMat(:,ib)==1,1,'first');
            if ~ismember(yUn,ytplc)
                ytplc = [ytplc yUn];
                ytlab = [ytlab Boundaries(ib)];
            end
        end
        set(gca,'ytick',fliplr(ytplc),'yticklabel',fliplr(ytlab))
    end
    if ist==1
%         ylabel('Unit #   |   z-scored FR, norm within each unit')
        xlabel('Time from stim onset')
    end
end %ist 


% colorbar('east')


% Save
print_eps_kp(hf1,fullfile(savedir,['PopOnset_' useFR '_' SortBy ]))



end





