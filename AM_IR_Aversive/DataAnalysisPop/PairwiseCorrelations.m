function PairwiseCorrelations(SUBJECT,SESSION)
%
%  allRastersSession(SUBJECT, SESSION )
%    Trying to check if neurons with similar tuning tend to have negative
%    correlations across trials. (after reading Yildiz 2016)
%    Might be better off defining similarity of tuning according to
%    spectral selectivity. Use response phase/type as a proxy? 
%
%  KP, 2019-09
%


%% Load files

fn = set_paths_directories(SUBJECT,SESSION);

filename = sprintf( '%s_sess-%s_Info'     ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_TrialData',SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,SUBJECT,SESSION); load(fullfile(fn.processed,SUBJECT,filename));

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Sort units by FR
[~, inspks] = sort(cellfun(@(x) sum(x>(TrialData.onset(1)/1000) & x<(TrialData.offset(1)/1000)), {Clusters.spikeTimes},'UniformOutput',true)); % baseline FR (during silence)
% [~, inspks] = sort(cellfun(@length, {Clusters.spikeTimes})); %overall nspikes
Clusters = Clusters(inspks);


%% Prepare figures

close all

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
rasterdotsize  = 6;
psthlinewidth  = 4;
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];
largerect = [1 scrsz(4)/2 scrsz(3)/3 scrsz(4)/2];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];


% Set colors
colors = [   0 200 150;...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;
colors = [ colors; ...
    [37  84 156]./255 ;...
    [19 125 124]./255 ;...
    [ 0   0   0]./255];

raster_colors = [0 0 0; 0.2 0.5 0.8];
patch_colors  = [0.9 0.9 0.9; 1 1 1];

Info.stim_ID_key{9} = 'Silence';
rmsrange = [];

% Get stimulus params
[dBSPL,LP] = theseSoundParams(TrialData);
if numel(dBSPL)>1 || numel(LP)>1
    keyboard
end


%%

Stimuli = 1:9;

TC_Corr   = nan(numel(Clusters),numel(Clusters));
NoiseCorr = nan(numel(Clusters),numel(Clusters));

for iUn1 = 1:numel(Clusters)
    for iUn2 = (iUn1+1):numel(Clusters)
        
        % Get spiketimes
        thisClu1 = Clusters(iUn1);
        spiketimes_1 = round(thisClu1.spikeTimes*1000 - spkshift)';
        
        thisClu2 = Clusters(iUn2);
        spiketimes_2 = round(thisClu2.spikeTimes*1000 - spkshift)';
        
        
        % Get FR stream
        bs_smth = 20;
        [Stream_FR_1,~,Stream_Spikes_1,~] = convertSpiketimesToFR(spiketimes_1,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        [Stream_FR_2,~,Stream_Spikes_2,~] = convertSpiketimesToFR(spiketimes_2,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        
%         [xc,lags]=xcorr(Stream_Spikes_1,Stream_Spikes_2,500);
        
        % Get all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials] = get_clean_trials(TrialData,[Info.artifact(thisClu1.maxChannel).trials; Info.artifact(thisClu2.maxChannel).trials],dBSPL,LP,1);
        allStim = unique(TrialData.trID(all_TDidx))';
        
        
        % Preallocate
        FR_mean = nan(numel(Stimuli),2);
        FR_sem  = nan(numel(Stimuli),2);
        FR_std  = nan(numel(Stimuli),2);
        
%         corr_withinTr  = nan(numel(Stimuli),1);
        corr_acrossTr  = nan(numel(Stimuli),1);
        
        for ist = 1:numel(Stimuli)
            
            stid = Stimuli(ist);
            
            Silence = false;
            if stid==9
                Silence = true;
            elseif ~ismember(stid,allStim)
                continue
            end
            
            ttp = min(Ntrials(Ntrials>0));
            
            if ist==1
                corr_withinTr  = nan(numel(Stimuli),ttp);
            end
            
            
            %% Collect trial indices and timestamps
            
            if ~Silence                                          % SOUND TRIALS
                
                if stid==3 || stid==6
                    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
                    % Find Pdc trials that follow same rate during ITI
                    TDidx = TDidx(TrialData(TDidx-1,:).trID ~= stid);
                    
                    TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
                    TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
                else
                    TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid );
                    TDidx_iti = [];
                end
                
                % Get timestamps of onsets and offsets
                clear t2 t3 Duration t_win
                t2 = TrialData.onset(TDidx);
                t3 = TrialData.offset(TDidx);
                Duration = mode(diff([t2 t3],1,2));
                
                % Add ITI trials (shortened to match duration)
                if ~isempty(TDidx_iti)
                    t2 = [t2; TrialData.onset(TDidx_iti)];
                    TDidx = [TDidx; TDidx_iti];
                end
                
            else                                               % SILENCE TRIALS
                
                clear t2 t3 Duration t_win
                SilPd = [TrialData.onset(1) TrialData.offset(1)];
                Duration = 1000;
                
                t2 = TrialData.onset(1) : Duration : (TrialData.offset(1)-mod(diff(SilPd),1000)-1000);
                TDidx = 1:length(t2);
                
            end
            
            
            %         kt = 1:length(t2);
            kt = sort(randperm(length(t2),ttp));
            t2     = t2(kt);
            TDidx  = TDidx(kt);
            
            t3 = t2 + Duration;
            
            if Silence && t3(end)>TrialData.offset(1)
                keyboard
            end
            
            
            %% Get data
            
            nSpk_tr   = nan( numel(TDidx), 2 );
            
            % Collect spikes/FR/rms for this stimulus/unit
            for it = 1:numel(t2)
                
%                 corr_withinTr(ist,it) = corr(Stream_Spikes_1(1,t2(it):t3(it))',Stream_Spikes_2(1,t2(it):t3(it))', 'type','Kendall');
                corr_withinTr(ist,it) = corr(Stream_FR_1(1,t2(it):t3(it))',Stream_FR_2(1,t2(it):t3(it))', 'type','Kendall');
                
                nSpk_tr(it,:) = [sum(Stream_Spikes_1(1,t2(it):t3(it))) sum(Stream_Spikes_2(1,t2(it):t3(it)))];
                
            end %it
                        
            r_acrossTr = corr(nSpk_tr,'type','Kendall');
            corr_acrossTr(ist,1) = r_acrossTr(1,2);
            
            nSpk_tr = nSpk_tr/Duration*1000;
            
            FR_mean(ist,:) = mean(nSpk_tr,1);
            FR_sem(ist,:)  = var(mean(nSpk_tr,1))./mean(nSpk_tr,1);
            FR_std(ist,:)  = std(nSpk_tr,1);
            
            
        end %ist
        
        
        SpTiming_Corr = mean(corr_withinTr,2,'omitnan');
        
        figure;
        plot([-0.5 0.5],[0 0],'k')
        hold on
        plot([0 0],[-0.5 0.5],'k')
        scatter(SpTiming_Corr,corr_acrossTr,5*sum(~isnan(corr_withinTr),2),colors,'filled')
        axis square
        xlabel('Spk timing corr')
        ylabel('Corr Nspk across trials')
        
        title(sprintf('Clus %i & %i',thisClu1.clusterID,thisClu2.clusterID))
        
        
        % Pair summary
        
        foo = corr(FR_mean(2:6,:));
        TC_Corr(iUn1,iUn2)   = foo(1,2);
        
        NoiseCorr(iUn1,iUn2) = mean(corr_acrossTr(2:6));
                
        
        % Plot
        
%         figure;
%         set(gcf,'Position',tallrect)
%         
%         % tuning curves
%         subplot(3,1,1)
%         plot(FR_mean(:,1),'c','LineWidth',2)
%         hold on
%         plot(FR_mean(:,2),'g','LineWidth',2)
%         xlim([0 ist+1])
%         set(gca,'xtick',1:numel(Stimuli),'xticklabel',Info.stim_ID_key(Stimuli))
%         title('tuning curves')
%         
%         % noise corr (across trials)
%         subplot(3,1,2)
%         plot([0 ist+1],[0 0],'k')
%         hold on
%         plot(corr_acrossTr,'k-o','LineWidth',2)
%         xlim([0 ist+1])
%         set(gca,'xtick',1:numel(Stimuli),'xticklabel',Info.stim_ID_key(Stimuli))
%         title('Noise Corr: FR across trials')
%         
%         % corr within each trial
%         subplot(3,1,3); 
%         plot([0 ist+1],[0 0],'k')
%         hold on
%         plotSpread(corr_withinTr')
%         boxplot(corr_withinTr')
%         xlim([0 ist+1])
%         set(gca,'xtick',1:numel(Stimuli),'xticklabel',Info.stim_ID_key(Stimuli))
%         title('time series correlation within each trial')
%         
%         suptitle(sprintf('Relationship of clus %i & %i',thisClu1.clusterID,thisClu2.clusterID))
        
    end % iUn2
end %iUn1


figure; 
plot(TC_Corr(:),NoiseCorr(:),'ok')
axis square

[r,p]=corr(TC_Corr(~isnan(TC_Corr(:))),NoiseCorr(~isnan(NoiseCorr(:))))


end