function AcceptMerge = checkUnits( Subject, Sessions, Channel, Clus, newSessLabel)
%  
%  called by setSplitUnits 
%  
%  Check some things before merging data from a unit that was
%  recorded across 2 sessions. Data will be concatenated at the
%  level of UnitData. 
%  
%  Unfinished alternative: merge processed data from sessions into a new
%  session.
%  
%  KP, 2018-12
% 

fn = set_paths_directories;

nClus = numel(Clus);
if nClus>2, keyboard, end

for iu = 1:nClus
    
    fprintf('Loading sess %s...\n',Sessions{iu})
    
    q=load(fullfile(fn.processed,Subject,sprintf( '%s_sess-%s_Info',     Subject,Sessions{iu})));
    eval(sprintf( 'Info%i = q.Info;', iu ))
    q=load(fullfile(fn.processed,Subject,sprintf( '%s_sess-%s_TrialData',Subject,Sessions{iu})));
    eval(sprintf( 'TrialData%i   = q.TrialData;' ,  iu ))
    eval(sprintf( 'SpoutStream%i = q.SpoutStream;', iu ))
    eval(sprintf( 'SoundStream%i = q.SoundStream;', iu ))
    q=load(fullfile(fn.processed,Subject,sprintf( '%s_sess-%s_Spikes'   ,Subject,Sessions{iu})));
    eval(sprintf( 'Spikes%i = q.Spikes;', iu ))
    
    clu_subtitle{iu} = sprintf('Sess %s: ch %i, clu %i ', Sessions{iu}, Channel, Clus(iu));
    
end

close all


%% Plot intepolated waveforms

for iu = 1:nClus
    
    eval( sprintf( 'spikes = Spikes%i.sorted(Channel);', iu ) )
    waveform = median(spikes.waveforms(spikes.assigns==Clus(iu),:),1);
    waveform = waveform./abs(min(waveform));
    m = 40;
    x = (1:length(waveform)) /Info1.fs*1000;
    q = linspace(min(x),max(x),length(x)*m);
    waveform_interp(iu,:) = interp1(x,waveform,q,'spline');
    
end

hf_w = figure;
plot(waveform_interp','LineWidth',3)
legend(clu_subtitle,'Location','northwest')
    


%% Plot autocorrelograms 

hf_a    = figure;
ISIvec = linspace(-100,100,100);

for iu = 1:nClus
    
    eval( sprintf( 'spikes = Spikes%i.sorted(Channel);', iu ) )
    ISIdist = diff( round(spikes.spiketimes(spikes.assigns==Clus(iu)).*1000) );
    ISIdist = [-1*ISIdist ISIdist];
    
    subplot(nClus,1,iu);
    hist(ISIdist,ISIvec)
    xlim([min(ISIvec)+20 max(ISIvec)-20])
    title(clu_subtitle{iu})
    
end

suptitle('Autocorrelograms')



%% (PCA plot (see AM stream folder) uses UMS program)



%% Plot rasters of responses to 4 & 32 Hz

these_stim = [3 6];

for ist = 1:numel(these_stim)
    
    hf_r(ist) = figure;
    
    for iu = 1:nClus
        
        eval( sprintf( 'TrialData = TrialData%i;', iu ))
        eval( sprintf( 'Info = Info%i;', iu ))
        
        % Get stimulus params
        [dBSPL,LP] = theseSoundParams(TrialData);
        if numel(dBSPL)>1 || numel(LP)>1
            keyboard
        end
        
        % Get spiketimes
        eval( sprintf( 'spikes = Spikes%i.sorted(Channel);', iu ) )
        spiketimes = round(spikes.spiketimes(spikes.assigns==Clus(iu)') * 1000);  %ms
        
        % Get trials
        all_TDidx = get_clean_trials(TrialData,Info.artifact(Channel).trials,dBSPL,LP);
        allStim = unique(TrialData.trID(all_TDidx))';
        
        stid = allStim(these_stim(ist));
        TDidx = all_TDidx( TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==0 );
        TDidx(TDidx==1) = [];
        
        TDidx_iti = [];
%         TDidx_iti = all_TDidx(TrialData.trID(all_TDidx)==stid & TrialData.ITIflag(all_TDidx)==1 & TrialData.Spout(all_TDidx)>0.95);
%         TDidx_iti = TDidx_iti(TrialData(TDidx_iti-1,:).trID>6);
        
        
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
        t3 = t2 + Duration;
        
        % Collect spikes/FR/rms for this stimulus/unit
        raster_x = [];
        raster_y = [];
        for it = 1:numel(TDidx)
            sp=[]; sp = spiketimes( spiketimes>=t2(it) ...
                & spiketimes<=t3(it) ) - t2(it) - 1;
            raster_x = [raster_x sp];
            raster_y = [raster_y it*ones(1,numel(sp))];
        end %it
        
        % Plot
        if numel(raster_x)>0
            
            subplot(nClus,1,iu);
            plot(raster_x, raster_y,...
                '.','MarkerSize',2)%,'Color',colors(stid,:))
            set(gca,'Color','none','ytick',it,'yticklabel',it,'ylim',[0 it+1])
            box off
            title(clu_subtitle{iu})
        end
        
    end %iu
    
    suptitle(['Raster plots, ' Info.stim_ID_key{these_stim(ist)}])
    
end %ist


%%
[AcceptMerge] = input('  1: Merge units.  0: Don''t merge.  ');


%% MERGE UNITS: CREATE NEW DATA FILES

if AcceptMerge
    
    newSessLabel
    
    
    length(SoundStream1)
    
    % Here, could make all files new & merged, even Spikes struct.
    % Be sure to check all programs for how they'd be impacted.
    
    
    
end



end