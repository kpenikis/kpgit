function PopRespSpeech
% 
% PopRespSpeech
% 
% 


%%% Find segments of Sound Stream with high cross correlation
%%% with threshold at 0.95, found ~52 instances of As You 
% figure;
% plot(Template_AsYou./max(Template_AsYou),'k','LineWidth',4)
% hold on
% 
% its   = 1:length(Template_AsYou);
% istep = 50;
% nlags = 100;
% 
% xc=0;
% while all(xc<0.95)
%     its = its+istep;
%     xc = crosscorr(Template_AsYou./max(Template_AsYou),SoundStream(its)./max(SoundStream(its)),nlags);
% end
% [~,im]=max(xc);
% SSseg = SoundStream((its(1)+im-nlags-1)+[1:length(Template_AsYou)]);
% 
% plot(SSseg./max(SSseg),'LineWidth',2)

close all

RERUN = 0;

bin_smooth = 20;
convwin     = 10;
convwin     = ones(1,convwin).*(1/convwin);

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
% spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load RepSpeech templates
k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
kfns = fieldnames(k);

Duration = max(structfun(@length,k));
Segs     = 1:numel(structfun(@length,k));
SegDurs  = structfun(@length,k);


savedir = fullfile(fn.figs,'PopulationSpeechSegments');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%%
if RERUN
    
    % Folded responses for repeated segments, from Pdc context only
    FR_vec    = nan(numel(UnitData),Duration,length(Segs));
    zFR_vec   = nan(numel(UnitData),Duration,length(Segs));
    
    for iUn = 1:numel(UnitData)
                
        % - - - -   Collect data   - - - - -
        
        get_trial_speech_RepSegs
        
        FR_vec(iUn,:,:)  = iu_FRvec;
        zFR_vec(iUn,:,:) = iu_zFRvec;
        % Save spiking data
        sp_trs(iUn,:,:,:) = thisRaster;
        
    end %iUn
    
    % Save MPH data
    save(fullfile(savedir,'RepSpeech'),'zFR_vec','FR_vec','RMS','sp_trs','-v7.3')
    
else
    
    
    
    %% Instead of re-running, load saved MPH data
    
    load(fullfile(savedir,'RepSpeech'))
    
    % Load an Info file too
%     filename = sprintf( '%s_sess-%s_Info',     UnitData(1).Subject,UnitData(1).Session); load(fullfile(fn.processed,UnitData(1).Subject,filename));
    
end



%% Plot results

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',12)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
twothirdsscreen  = [1 scrsz(4) scrsz(3)/3*2 scrsz(4)];
matchedsize      = [1 scrsz(4)/2 scrsz(3)/5*2 scrsz(4)/2];

spwidth = [1 3 6 9 10];


%%
useFR = 'log';

SortBy =  'zThLatAbv'; 'latency'; 'signsum'; 'mean'; 'range'; 
thresh = 0.2;
inWin  = 0.25;

%~~~~~~~~~~~~
%   Fig 1      = PSTH
%~~~~~~~~~~~~
hf1 = figure;
set(gcf,'Position',matchedsize)

%~~~~~~~~~~~~
%   Fig 2      = rasters
%~~~~~~~~~~~~
nTrs = 5;
for itr = 1:nTrs
    hf2(itr)=figure;
    set(gcf,'Position',matchedsize)
    hold on
end

for ist=1:size(FR_vec,3)
    
    switch useFR
        case 'z'
            ThisData    = zFR_vec(:,1:SegDurs(ist),ist);
            Boundaries  = [-1 0 0.25 0.5 1 2];
            SpData      = sp_trs(:,1:SegDurs(ist),ist,:);
        case 'log'
            ThisData    = FR_vec(:,1:SegDurs(ist),ist);
            Boundaries  = [-1 round(10.^[0.5 1 1.25 1.5]) ];
            SpData      = sp_trs(:,1:SegDurs(ist),ist,:);
    end
    
    % Sort here
    switch SortBy
        case 'mean'
            [~,i_sorted] = sort(mean(ThisData,2));
        case 'range'
            [~,i_sorted] = sort(range(ThisData,2));
        case 'signsum'
            [~,i_sorted] = sort(sum(sign(ThisData),2)) ; 
        case 'zThLatAbv'
            sortdata = [];
            for iUn=1:size(ThisData,1)
                clear uThZ uLat uAbv
                if all(isnan(ThisData(iUn,:,1)))
                    sortdata = [sortdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = Boundaries(find(max(ThisData(iUn,:,1),[],2) > Boundaries,1,'last'));
                    %find time exceeded it
                    uLat = find(ThisData(iUn,:,1) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(ThisData(iUn,:,1) > uThZ);
                    
                    % Save unit data
                    sortdata = [sortdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( sortdata, [1 -2 -3]);
            
        case 'latency'
            zTh = [-0.5 0.5];
            sortdata = [];
            for iUn=1:size(ThisData,1)
                clear uThZ uLat uAbv
                if all(isnan(ThisData(iUn,:,1)))
                    sortdata = [sortdata; nan nan nan ];
                else
                    %find highest thresh exceeded
                    uThZ = zTh(find(max(ThisData(iUn,:,1),[],2) > zTh,1,'last'));
                    %find time exceeded it
                    uLat = find(ThisData(iUn,:,1) > uThZ,1,'first');
                    %find ms above
                    uAbv = sum(ThisData(iUn,:,1) > uThZ);
                    
                    % Save unit data
                    sortdata = [sortdata; uThZ uLat uAbv ];
                end
            end
            [~,i_sorted] = sortrows( sortdata, [1 -2]);
    end
    
    %~~~~~~~~~~~~
    %   PSTH
    %~~~~~~~~~~~~
    figure(hf1); hold on
    
    % Plot stim RMS
    subplot(4,9,spwidth(ist):(spwidth(ist+1)-1));
    
    plot(RMS(ist,:),'k','LineWidth',2)
    xlim([0 SegDurs(ist)])
    ylim([0 0.04])
    set(gca,'ytick',[],'xtick',[],'Color','none')
    title(kfns{ist})
    box off
    
    % Plot responses
%     subplot(4,size(FR_vec,3),ist+size(FR_vec,3).*[1 2 3]);
    spvec = (spwidth(ist):(spwidth(ist+1)-1))'+9.*[1 2 3];
    subplot(4,9,sort(spvec(:))');
    
    plotdata = ThisData(i_sorted,:,1) ;
    
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
            cmocean('-gray')
    end
    
    xlim([0 SegDurs(ist)])
    ylim([0.5 ndp+0.5])
    set(gca,'ytick',[],'xtick',[0 SegDurs(ist)],'tickdir','out','ticklength',[0.02 0.02],'Color','none')
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
    
    
    
    %~~~~~~~~~~~~
    %  Rasters
    %~~~~~~~~~~~~
        
    for itr = 1:nTrs
        
        plotspks = SpData(i_sorted,:,:,itr);
        
        figure(hf2(itr));
        hold on
        
        % Plot stim RMS
        subplot(4,9,spwidth(ist):(spwidth(ist+1)-1));
        
%         plot(RMS(ist,:),'k','LineWidth',2)
%         xlim([0 SegDurs(ist)])
%         ylim([0 0.04])
%         set(gca,'ytick',[],'xtick',[],'Color','none')
%         title(kfns{ist})
%         box off
        
        % OR plot PSTH
        plot(sum(plotspks(1:ndp,:),1),'k')
        smFR = smoothFR(sum(plotspks(1:ndp,:),1)./ndp,20);
        hold on
%         plot(smFR,'b')
        foo = conv(sum(plotspks(1:ndp,:),1),convwin);
        ps_conv = foo(floor(numel(convwin)/2)+(0:size(plotspks,2)-1));
        plot(1:length(ps_conv),ps_conv,'m','LineWidth',2)
        
        xlim([0 SegDurs(ist)])
        set(gca,'xtick',[],'Color','none')
        title(kfns{ist})
        
        % Plot spikes
        subplot(4,9,sort(spvec(:))');
        hold on
        
        % raster
        [x,y] = find(plotspks(1:ndp,:)');
        plot([x x]',y'+[-0.5; 0.5],'-k')
        set(gca,'ydir','reverse')
        
        % Add markers to label NS cells
        %         flagNS = UnitInfo.TroughPeak(i_sorted(1:ndp))<0.5;
        %         hold on
        %         plot(0,find(flagNS),'.','Color',[0.01 0.57 0.44])
        %         plot(2*ceil(1000/AMrates(ir)),find(flagNS),'.','Color',[0.01 0.57 0.44])
        
        % Finish plot
        xlim([0 SegDurs(ist)])
        ylim([0.5 ndp+0.5])
        set(gca,'tickdir','out','ticklength',[0.02 0.02],'Color','none')
        set(gca,'xtick',[0 SegDurs(ist)],'ytick',[])
        %     [r,~]=find(cumsum(sortdata(i_sorted(1:ndp),1)==Boundaries)==1);
        %     set(gca,'ytick',r,'yticklabel',Boundaries)
        
        box off
        axis fill
        
        set(gcf,'PaperOrientation','landscape')
        
        
    end %itr
    
    
end %ist 

% colorbar('east')


% Save PSTHs
savename = ['RepSpeech_' useFR '_' SortBy ];
print_eps_kp(hf1,fullfile(savedir,savename))
set(hf1,'PaperOrientation','landscape')
print(hf1,'-dpdf','-r500','-fillpage', fullfile(savedir,savename))


% Save rasters 
for itr = 1:nTrs
    savename = sprintf('%s_Tr%i',savename,itr);
    print_eps_kp(hf2(itr),fullfile(savedir,'Trials',savename))
%     print(hf2(itr),'-dpdf','-r500','-fillpage', fullfile(savedir,'Trials',savename) )
end



end



