function ap_unmod_raster(select_subject,select_session,select_channel,Clus)
%
%  ap_unmod_raster(subject,session,channel,clu)
%    
%
%  KP, 2017-09
%


close all

%!!!!!!!!!!!!!!!!!
SUonly   =  1;
%!!!!!!!!!!!!!!!!!
FRcutoff =  3;%Hz 
%!!!!!!!!!!!!!!!!!
minTrs   =  10;
%!!!!!!!!!!!!!!!!!

%%
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');
figsize1 = [1 scrsz(4)*3/4 scrsz(3)*2/3 scrsz(4)*3/4];
figsize1 = 1+0.7.*[1 scrsz(4)/5 scrsz(3)/6 scrsz(4)/5];


colors = [ 84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0;...
           255 205  60 ]./255;
colors = [colors; 0.7.*bone(4)];


%%
unType = {'' 'SU' 'MU'};

N = 0;

histbinsize = 20;
anbinsize   = 50;
smthbinsize = 50;


% Make empty data table
Zdata = table;
Zdata.Subject  = ' ';
Zdata.Session  = ' ';
Zdata.ch       = 0;
Zdata.clu      = 0;
Zdata.IRblock  = 0;
Zdata.stimpars = [0 0];
Zdata.zFR_pred = 0;
Zdata.zFR_obs  = 0;
% Zdata.preUNM   = 0;
% Zdata.postUNM  = 0;



%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    subjects = {select_subject};
else
    subjects = {'WWWf_244303' 'WWWr_244300'};
end

for subj = 1:numel(subjects)

    subject = subjects{subj};
    
    switch subject
        case 'WWWf_244303'
            subjcol = 'k';
        case 'WWWr_244300'
            subjcol = 'k';
    end
     

%%  SESSIONS

% Get list of sessions to check for sorted data

fn = set_paths_directories(subject);

if nargin>1 && exist('select_session','var')
    Sessions = {select_session};
else
    
    SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
    Sessions = [];
    for ifn = 1:numel(SpkFns)
        if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
            Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
        end
    end
    
end
% Sessions = flipud(Sessions);

% Step through each session
for sess = Sessions'
    
session = char(sess);

%%
% Load data files
fn = set_paths_directories(subject,session);
fprintf('Session %s:  loading data...\n',session)
filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
filename = sprintf( '%s_sess-%s_SoundData',subject,session); load(fullfile(fn.processed,subject,filename));
load(fullfile(fn.stim,'IRsequences.mat'))


%%
% GET STIM INFO

AMrates = [2 4 8 16 32 64];

TotalDur_s = size(SoundData,2)/Info.fs_sound;


% % Get unique dBSPLs
% dBSPL = unique(SoundData(4,:));
% rm_i=[];
% for ii = 1:numel(dBSPL)
%     if (numel(find(SoundData(4,:)==dBSPL(ii)))/Info.fs_sound) < 60
%         rm_i = [rm_i ii];
%     end
% end
% dBSPL(rm_i) = [];
% 
% % Get unique noisebands (based on LP)
% LP = unique(SoundData(6,:));
% rm_i=[];
% for ii = 1:numel(LP)
%     if (numel(find(SoundData(6,:)==LP(ii)))/Info.fs_sound) < 60
%         rm_i = [rm_i ii];
%     end
% end
% LP(rm_i) = [];
% 
% % Get unique AM depths
% AMdepth = unique(SoundData(3,:));
% rm_i=[];
% for ii = 1:numel(AMdepth)
%     if (numel(find(SoundData(3,:)==AMdepth(ii)))/Info.fs_sound) < 60
%         rm_i = [rm_i ii];
%     end
% end
% AMdepth(rm_i) = [];



%% STEP THROUGH EACH CHANNEL

if nargin>2 && exist('select_channel','var')
    Channels = select_channel;
else
    Channels =  [1:7 9:16];
end

for channel = Channels
        
    % Artifact for this channel
    ArtifactFlag = Info.artifact(channel).SDsamples;
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if all(spikes.labels(:,2)==1)
            continue
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
            continue
        else
            Clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = Clus'
        
        % !! Only SU for now !!
        if SUonly && spikes.labels(spikes.labels(:,1)==clu,2) ~= 2
            continue
        end
        
        try
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        
        % Skip units with overall FR in session below a cutoff (3 hz)
        if numel(spiketimes) < round(FRcutoff*length(SoundData)/Info.fs_sound)
            continue
        end
        
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~
        % Convert FR to z-score
        %~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Make constant stream of 
        Stream_Spks = zeros(1,1000*ceil((size(SoundData,2)/Info.fs_sound)));
        Stream_Spks(spiketimes) = 1;
        
        %either with standard 20 ms bin
%         Stream_FRbin = 1000*(binspikecounts(Stream_Spks,histbinsize)/histbinsize);
%         Stream_FRbin(isinf(Stream_FRbin)) = nan;
%         foo = repmat(Stream_FRbin,histbinsize,1);
%         Stream_FR = reshape(foo,1,histbinsize*length(Stream_FRbin));
%         Stream_FR = Stream_FR(1:ceil(length(SoundData)/Info.fs_sound*1000));
        
        %or with sliding 50 ms boxcar
        Stream_FRsmooth = smoothFR(Stream_Spks,smthbinsize);
        Stream_FRsmooth = Stream_FRsmooth(1:ceil(length(SoundData)/Info.fs_sound*1000));
        
        % Determine the time of the beginning of actual recording session
        % max of first spiketime or 5 seconds before unmodulated sound came
        % on
        sampStart = find(diff(SoundData(8,:))==-1);
        msStart   = max( spiketimes(1), round(sampStart(1)/Info.fs_sound*1000)-5000 );
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Convert FR to z-score
        
        REFERENCE = 'unmod';'stream';
        markRED = 0;
        
        switch REFERENCE
            
            case 'stream'
                Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
            
            case 'unmod'
                
                if ~any(strcmp(Zdata.Properties.VariableNames,'preUNM'))
                    Zdata.preUNM   = 0;
                    Zdata.postUNM  = 0;
                end
                
                % Get samples of unmodulated sound and of silence
                unmod_samps  = find(SoundData(8,:)==11);
                
                % Split unmod portion into before and after
                try
                    if ~isempty(find(diff(unmod_samps)>1))
                        unmod1_ms = round( [unmod_samps(1)   unmod_samps(0+find(diff(unmod_samps)>1))] /Info.fs_sound*1000 );
                        unmod2_ms = round( [unmod_samps(1+find(diff(unmod_samps)>1)) unmod_samps(end)] /Info.fs_sound*1000 );
                        
%                         diff(unmod1_ms)/1000
%                         diff(unmod2_ms)/1000
                        
                        meanFR = mean([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        stdFR = std([Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)) Stream_FRsmooth(unmod2_ms(1):unmod2_ms(2))]);
                        
                    else
                        
                        warning(' !! Can''t make a raster for this unit, because there wasn''t unmodulated sound at the end of the recording.')
                        keyboard
                        
                        unmod1_ms = round( [unmod_samps(1) unmod_samps(end)]/Info.fs_sound*1000 );
                        unmod2_ms = [0 0];
                        
                        meanFR = mean(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        stdFR = std(Stream_FRsmooth(unmod1_ms(1):unmod1_ms(2)));
                        
                        fprintf(' !! not much unmod data for sess %s ch%i clu%i \n',session,channel,clu)
%                         markRED = 1;
                    end
                    
                catch
                    keyboard
                end
                
                Stream_zscore = (Stream_FRsmooth(msStart:end) - meanFR) / stdFR;
                Stream_zscore = [nan(1,msStart-1) Stream_zscore];
                
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
                
        % Get unique dBSPLs, noisebands (based on LP), and AM depths
        dBSPL = unique(SoundData(4,unmod_samps));
        LP = unique(SoundData(6,unmod_samps));
        AMdepth = unique(SoundData(3,unmod_samps));
        
        if any([numel(dBSPL) numel(LP) numel(AMdepth)]>1)
            keyboard
        end
        
        
        %%
        
        % Get spiking data from arbitrary chunks of unmod data
        
        % Define chunks for sound at end
        chunks = [unmod2_ms(1):diff(unmod1_ms):unmod2_ms(2)];
        
        % Create empty arrays for data
        raster_x  = [];
        raster_y  = [];
        stim      = nan(numel(chunks),ceil(diff(unmod1_ms)/1000*Info.fs_sound));
        PSTH      = nan(numel(chunks),diff(unmod1_ms));
        
        % First get spikes for sound preceding AM stream
        chunk_samps = [floor(unmod1_ms(1)/1000*Info.fs_sound) ceil((unmod1_ms(2)-1)/1000*Info.fs_sound)];
        raster_x = [raster_x spiketimes( spiketimes>=unmod1_ms(1) & spiketimes<unmod1_ms(2) ) - unmod1_ms(1) ];
        raster_y = [raster_y ones(1,sum(spiketimes>=unmod1_ms(1) & spiketimes<unmod1_ms(2))) ];
        stim(1,1:diff(chunk_samps)+1) = envelope( SoundData(2,chunk_samps(1):chunk_samps(2)) ,50,'rms');
        PSTH(1,:) = Stream_FRsmooth(unmod1_ms(1):(unmod1_ms(2)-1));
        
        % Then get spikes for sound at end        
        for it = 2:numel(chunks)
            chunk_samps = [floor(chunks(it-1)/1000*Info.fs_sound) ceil((chunks(it)-1)/1000*Info.fs_sound)];
            sp=[];  sp  = spiketimes( spiketimes>=chunks(it-1) & spiketimes<chunks(it) );
            raster_x    = [raster_x sp-chunks(it-1) ];
            raster_y    = [raster_y it*ones(1,numel(sp))];
            PSTH(it,:)   = Stream_FRsmooth(chunks(it-1):(chunks(it)-1));
            stim(it,1:diff(chunk_samps)+1)  = envelope( SoundData(2,chunk_samps(1):chunk_samps(2)) ,50,'rms');
        end
        
        if any(sum(isnan(PSTH),1))
            keyboard
        end
        
        % RASTER/PSTH
        
        figure;
        set(gcf,'Position',figsize1)
        
        % STIMULUS
        subplot(5,1,1)
        plot(stim','Color',[0.5 0.5 0.5])
        xlim([0 size(stim,2)])
        set(gca,'xtick',[],'ytick',[])
%         ylabel('stim rms')
        set(gca,'FontSize',12)
        
        %         titlestr1 = sprintf('unmodulated noise (reference stimulus)\n%s | %s | ch %i, clu %i (%s)\n%i dBSPL | %i-%i Hz | depth=%i',...
        %             subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},dBSPL,SoundData(5,unmod1_ms(1)),LP,AMdepth);
        %         title(titlestr1)

        % RASTER
        subplot(5,1,2:3); hold on
        scatter(  raster_x  ,  raster_y,  3, 'o',...
            'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none')
        ylim([0 it+1])
        xlim([0 diff(unmod1_ms)])
        set(gca,'xtick',[],'ytick',it+1,'FontSize',12)
        ylabel('chunks')        
        
        % PSTH
        subplot(5,1,4:5); hold on
        plot(mean(PSTH,1,'omitnan'),'Color',[0.5 0.5 0.5],'LineWidth',3)
        ylim([0 120])
        xlim([0 diff(unmod1_ms)])
        xlabel('time (ms)')
%         ylabel('spikes/s')
        set(gca,'FontSize',12,'xtick',[],'ytick',120)
        
        
        % Save figure
        
        savedir = fullfile(fn.rasterplots,['ch' num2str(channel)]);
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        savename = sprintf('%s_%s_ch%i_clu%i_%s_UNMOD_%idB_LP%ihz',...
            subject,session,channel,clu,unType{spikes.labels(spikes.labels(:,1)==clu,2)},dBSPL,LP);
        print_eps_kp(gcf,fullfile(savedir,savename),1)
        
%         set(gcf,'PaperOrientation','landscape');
%         print(gcf,'-dpdf',fullfile(savedir,savename),'-bestfit')
%         print(gcf,fullfile(savedir,savename),'-depsc','-tiff')
        
        
        
    end %clu
end %channel

end % sessions

end % subject


end %function

%}


