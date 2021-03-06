function allSession_zScore(SUBJECT,SESSION)
%
%  allRastersSession(SUBJECT, SESSION )
%    Plots a raster and psth for each stimulus, for all SU from one session.
%    Uses Clusters struct, not UnitData table. 
%
%  KP, 2019-01, updated 2019-02
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
set(0,'DefaultAxesFontSize',12)
psthlinewidth  = 2;
rng('shuffle')

scrsz = get(0,'ScreenSize');   %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

figwidthscales = [1.5 1 1 1 1 1 1.937 1.937 1];

Info.stim_ID_key{9} = 'Silence';
rmsrange = [];


%%
Stimuli = 1:9;
for ist = Stimuli
    
    stid = Stimuli(ist);
    
    add_y    = 0;%zeros(1,7);
    N_un     = 0;%zeros(1,7);
    ytickstr = cell(1,numel(Clusters));
    zDataALL = [];
    
    %%
    % Set up figure
    hf(ist) = figure;
    set(gcf,'Position',tallrect.* [1 1 figwidthscales(stid) 1])
    hold on
    hs(ist,1) = subplot(9,1,1);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[],'ytick',[])
    hs(ist,2) = subplot(9,1,2:9);
    set(gca,'xlim',[0 figwidthscales(stid)*1000],'xtick',[0 (figwidthscales(stid)*1000)],'ytick',[],...
        'xticklabel',[0 (figwidthscales(stid)*1000)])
    
    
    for iUn = 1:numel(Clusters)
        
        % FOR NOW, ADD PLACEHOLDER FOR ARTIFACT TRIALS FOR SYNAPSE/KS DATA
        if ~isfield(Info,'artifact')
            keyboard
            Info.artifact(64).trials = [];
        end        
        
        % Get spiketimes
        thisClu = Clusters(iUn);
        maxChan = thisClu.maxChannel;
        %originally: seconds, not rounded
        spiketimes = unique(round(thisClu.spikeTimes*1000-spkshift)');
        
        
        % Get stimulus params
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        
        % Convert FR to z-score
        bs_smth = 20;
        [Stream_FRsmooth,Stream_zscore,Stream_Spikes,ymaxval] = convertSpiketimesToFR(spiketimes,...
            length(SpoutStream),TrialData.onset(1),TrialData.offset(1),10,bs_smth,'silence');
        
        meanFR = 1000* sum(Stream_Spikes(TrialData.onset(1):end)) / length(Stream_Spikes(TrialData.onset(1):end));
        
        
        % Get all stimuli presented with these parameters, given a
        % sufficient number of trials without diruptive artifact
        % while the animal was drinking
        [all_TDidx,Ntrials] = get_clean_trials(TrialData,Info.artifact(maxChan).trials,dBSPL,LP);
        allStim = unique(TrialData.trID(all_TDidx))';
        
        Silence = false;
        if stid==9
            Silence = true;
        elseif ~ismember(stid,allStim)
            continue
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
        
        
        kt = 1:length(t2);
        t2     = t2(kt);
        TDidx  = TDidx(kt);
        
        t3 = t2 + Duration;
        
        if Silence && t3(end)>TrialData.offset(1)
            keyboard
        end
        
        
        %%
        % Preallocate
        zData  = nan( numel(TDidx), Duration+1 );
        stim   = nan( numel(TDidx), Duration+1 );
        
        % Collect spikes/FR/rms for this stimulus/unit
        for it = 1:numel(TDidx)
            
            stim(it,:) = ...
                SoundStream(1, t2(it) : t3(it) )...
                ./ max(SoundStream(1, t2(it) : t3(it) ));
            
            zData(it,:) = ...
                Stream_zscore(1, t2(it) : t3(it) );
            
        end %it
        
        if Silence
            stim = zeros(size(stim));
        end
        
        
        %% Add to plots
        
        figure(hf(ist)); hold on
        
        % Stimulus RMS
        subplot(hs(ist,1)); hold on
        plot(0:Duration, mean(stim,1,'omitnan'),...
            'LineWidth',psthlinewidth,'Color','k')
        set(gca,'Color','none','xtick',[],'ytick',[])
        
        % Save average zScore activity for this unit
%         subplot(hs(ist,2)); hold on
%         image(0:Duration, N_un*ones(1,size(zData,2)), mean(zData,1,'omitnan') )
        zDataALL = [zDataALL; mean(zData,1,'omitnan')];
        
        add_y = [add_y add_y(end) + max(kt)];
        N_un = N_un + 1;
        ytickstr{N_un} = num2str(thisClu.clusterID);
        
    end % iUn
    
    if add_y == 0
        continue
    end
    
    subplot(hs(ist,2)); hold on
    
    imagesc(zDataALL)
    
    cmocean('balance','pivot',0)
    caxis([-1 3])
    
    box off    
    set(gca,'Color','none','xtick',[],...
        'tickdir','out','ticklength',0.01.*[1 1])
    
    set(gca,'ytick',1:N_un,'yticklabel',ytickstr,'ylim',[0.5 N_un+0.5],...
        'xtick',unique([0:500:Duration Duration]),'xticklabel',unique([0:500:Duration Duration]))
    xlabel('Time (ms)')
%     ylabel([num2str(mode(diff(add_y))) ' trials each'])
    
    subplot(hs(ist,1));
    rmsrange = [rmsrange; get(gca,'ylim')];
    set(gca,'ylim',mode(rmsrange,1))
    
        
    suptitle(sprintf('%s   |   %s %s   |   %i units',...
        Info.stim_ID_key{stid}, SUBJECT,SESSION,N_un ))
    
    
    
    %% Save figure
    
    savedir = fullfile(fn.processed,'Rasters_allSess',sprintf('%s_%s',SUBJECT,SESSION));
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    savename = sprintf('%s_%s_zSc_%s',SUBJECT,SESSION,Info.stim_ID_key{ist});
%     print_eps_kp(hf(ist),fullfile(savedir,savename))
    set(gcf,'PaperOrientation','landscape')
    print(hf(ist),fullfile(savedir,savename),'-dpdf','-bestfit')
    
    
end %ist


% Print figure for legend
hcb=figure;
imagesc(zDataALL)
cmocean('balance','pivot',0)
caxis([-1 3])
colorbar
set(gca,'ytick',1:N_un,'yticklabel',ytickstr,'ylim',[0.5 N_un+0.5])

savename = sprintf('%s_%s_zSc_Colorbar',SUBJECT,SESSION,Info.stim_ID_key{ist});
% print_eps_kp(hcb,fullfile(savedir,savename))
set(gcf,'PaperOrientation','landscape')
print(hf(ist),fullfile(savedir,savename),'-dpdf','-bestfit')



end