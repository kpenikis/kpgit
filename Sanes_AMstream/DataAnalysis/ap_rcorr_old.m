function ap_rcorr(subject)
%
%  ap_rcorr(subject)
%    
%
%  KP, 2017-08
%

global mspad

% Add some time before and after block, to pad the convolution
mspad = 500;



%%  SESSIONS

% Get list of sessions to check for sorted data

fn = set_paths_directories(subject);
SpkFns = dir(fullfile(fn.processed,subject,'*_Spikes.mat'));
    
Sessions = [];
Sessions2 = [];
for ifn = 1:numel(SpkFns)
    splitStr = regexp(SpkFns(ifn).name,'_','split');
    splitStr2 = regexp(splitStr{3},'-','split');
    if length(splitStr2{2})==2
        Sessions = [Sessions; splitStr2{2}];
    end
%     if length(char(extractBetween(SpkFns(ifn).name,'sess-','_Spikes')))==2
%         Sessions2 = [Sessions2; extractBetween(SpkFns(ifn).name,'sess-','_Spikes')];
%     end
end


% Set some plotting params
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');  % [ left bottom width height ]
figsize = [1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)/2];

% Set save directory
savedir = fullfile(fn.processed,subject,'^an_plots','RCorr');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


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


% Get unique dBSPLs
dBSPL = unique(SoundData(4,:));
rm_i=[];
for ii = 1:numel(dBSPL)
    if (numel(find(SoundData(4,:)==dBSPL(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
dBSPL(rm_i) = [];

% Get unique noisebands (based on LP)
LP = unique(SoundData(6,:));
rm_i=[];
for ii = 1:numel(LP)
    if (numel(find(SoundData(6,:)==LP(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
LP(rm_i) = [];

% Get unique AM depths
AMdepth = unique(SoundData(3,:));
rm_i=[];
for ii = 1:numel(AMdepth)
    if (numel(find(SoundData(3,:)==AMdepth(ii)))/Info.fs_sound) < 60
        rm_i = [rm_i ii];
    end
end
AMdepth(rm_i) = [];




%% GET SPIKE TIMES

% Step through all channels if not specified
if nargin<3 && ~exist('channels','var')
    channels = [1:7 9:16];
end


%% STEP THROUGH EACH CHANNEL
for channel = channels
        
    % Artifact for this channel
    ArtifactFlag = Info.artifact(channel).SDsamples;
    
    % Find clus to plot
    spikes = Spikes.sorted(channel);
    if nargin<4
        if all(spikes.labels(:,2)==1)
%             disp(' SESSION MAY NOT BE MANUALLY SORTED YET')
            continue
        end
        if ~any(spikes.labels(:,2)==2 | spikes.labels(:,2)==3)
%             disp(' !! no valid clus for this channel')
            continue
        else
            clus = spikes.labels(spikes.labels(:,2)==2 |spikes.labels(:,2)==3, 1);
        end
    end
    
    %% STEP THROUGH EACH CLU
    
    for clu = clus'
        
        % !! Only SU for now !!
        if spikes.labels(spikes.labels(:,1)==clu,2) ~= 2
            continue
        end
        
        try
        spiketimes = round(spikes.spiketimes(spikes.assigns==clu') * 1000);  %ms
        catch
            keyboard
        end
        % spiketrials = spikes.trials(unit_in);
        
        if isempty(spiketimes)
            error('no spike events found for this clu')
            %also skip clus with few events 
        elseif numel(spiketimes) < round(3*length(SoundData)/Info.fs_sound)
            continue
        end
        
        
        
        %%        
        
        % Step through each combo of dBSPL, LP filter, AMdepth
        for spl = dBSPL
            for lpn = LP
                for amd = AMdepth
                    
                    % Get this unit/stimulus combo's N clean blocks (trials) 
                    [blocksN,minDur] = get_N_clean_blocks(SoundData,Info,ArtifactFlag,spl,lpn,amd);
                    
                    if any(blocksN<10), continue, end
                    
                    fprintf('calculating RCORR vals for ch %i clu %i...\n',channel,clu)
                    
                    %~~~~~~~~
                    sws = [1 2 4];% 8 16 32 64 128 256];
                    nTemps = 5;
                    %~~~~~~~~
                    
                    percentCorrectMat = nan(numel(sws),2,nTemps);
                    best_PC   = 0;
                    best_sw   = 0;
                    best_data = [];
                    
                    for isw = 1:numel(sws)
                        
                        % Create gaussian smoothing window for convolving spiketrains
                        sw = sws(isw);
                        [~,GW,~] = makeSmoothGauss_mlc(sw,minDur);
                        
                        % Determine the half length (in samples) of the gaussian window. Use later
                        % for shifting convolved signal to align with correct time points.
                        lGausHalf = (length(GW)-1)/2;
                        
                        % Create empty matrix to track the template trials used
                        ntUsed = zeros(1,length(blocksN));
                        
                        percentAssignment = nan(10,10,nTemps);
                        
                        for iTemp = 1:nTemps
                            
                            [T,ntUsed] = getTemplates(blocksN,minDur,ntUsed,spiketimes,...
                                spl,lpn,amd,ArtifactFlag,Info.fs_sound,SoundData,GW);
                            
                            
                            % Preallocations
                            finalMat = [];
                            count = 0;
                            
                            % Go through each block type
                            for thisblock = 1:numel(Info.blockKey)
                                
                                %if marked unmodulated or silent, skip
                                if thisblock>=11, continue, end
                                
                                % Get onsets for this block type
                                [~,~,bkStart_ms,bkStop_ms]  =  get_blockOnsets( SoundData,...
                                    thisblock,spl,lpn,amd,ArtifactFlag,Info.fs_sound);
                                
                                % Go through each trial
                                for it = 1:numel(bkStart_ms)
                                    
                                    % Skip the current template trial
                                    if it==ntUsed(end,thisblock)
                                        continue
                                    end
                                    
                                    % Get spiketimes
                                    sp=[];
                                    sp = spiketimes( spiketimes>=(bkStart_ms(it)-mspad) & spiketimes<(bkStop_ms(it)+mspad) ) - bkStart_ms(it) + 1 + mspad;
                                    
                                    % Convert from spiketimes to binary vector
                                    sp01 = zeros(1,bkStop_ms(it)-bkStart_ms(it)+2*mspad);
                                    sp01(sp) = 1;
                                    
                                    % Convolve with Gaussian window
                                    sp_conv = conv(sp01,GW,'full');
                                    
                                    % Correct for the time shift introduced by the convolution
                                    sp_conv = sp_conv(lGausHalf+1:end-lGausHalf);
                                    
                                    % Check by plotting convolved signal against spiketimes
                                    %     figure; hold on
                                    %     h = plot(sp_conv,'r-');
                                    %     y = zeros(length(sp),1);
                                    %     h2 = plot(sp,y,'r.');
                                    
                                    % Remove extra time added for padding AND cut block to min block dur
                                    S = sp_conv(mspad+(1:minDur));
                                    
                                    % Check by plotting again
                                    %     figure; hold on
                                    %     h = plot(sp_conv,'b-');
                                    %     y = zeros(length(sp),1);
                                    %     h2 = plot(sp-mspad,y,'b.');
                                    %     xlim([1 minDur])
                                    
                                    
                                    % Calculate RCORR values and make block assignment
                                    
                                    [blockAssignment,maxR] = calcR(T,S);
                                    count = count+1;
                                    finalMat(count,:) = [thisblock,it,blockAssignment,maxR];
                                    
                                    
                                end %it
                                
                                
                            end %thisblock
                            
                            
                            % Calculate overall percent correct
                            correct = find(finalMat(:,1) == finalMat(:,3));
                            percentCorrect = (numel(correct)/length(finalMat))*100;
                            percentCorrectMat(isw,:,iTemp) = [sw,percentCorrect];
                            
                            
                            % Collect data to visualize full results
                            for thisblock = 1:max(finalMat(:,1))
                                percentAssignment(thisblock,:,iTemp) = hist(finalMat(finalMat(:,1) == thisblock,3),1:10)./sum(finalMat(:,1) == thisblock);
                            end
                            
                            
                        end %iTemp
                        
                        
                        % If this sw has a higher overall PC than current,
                        % save the data to plot later in MTF
                        if mean(percentCorrectMat(isw,2,:),3)>best_PC
                            best_sw   = sw;
                            best_PC   = mean(percentCorrectMat(isw,2,:),3);
                            best_data = diag(mean(percentAssignment,3)).*100;
                        end
                            
                        
                        
                        % Plot full results, averaged across the random
                        % template trials chosen
                        hf=figure; hold off                        
                        imagesc(mean(percentAssignment,3)')
                        axis square
                        set(gca,'xtick',1:10,'ytick',1:10)
                        titlestr = sprintf('%s | %s | ch %i, clu %i\n%i dBSPL | %i-%i Hz\nRCorr results, with filter sigma %i ms',...
                            subject,session,channel,clu,spl,SoundData(5,bkStart_ms(1)),lpn,sw);
                        title(titlestr)
                        xlabel('True block')
                        ylabel('Assigned block')
                        colorbar
                        caxis([0 0.6])
                        
                        
                        % Save figure
                        savename = sprintf('%s_%s_ch%i_clu%i_%idB_LP%ihz_RCorr_sw%i',...
                            subject,session,channel,clu,spl,lpn,sw);
                        set(hf,'PaperOrientation','landscape');
                        print(hf,'-dpdf',fullfile(savedir,savename),'-bestfit')
                        
                        
                    end %isw
                    
                    
                    % And plot the MTF/tuning function of the best sw
                    hf1=figure; hold on
                    set(hf1,'Position',figsize)
                    plot(1:6,best_data(1:6),'k.--','MarkerSize',30,'LineWidth',1)
                    plot(7:10,best_data(7:10),'k.--','MarkerSize',30,'LineWidth',1)
                    set(gca,'xlim',[0 11],'ylim',[0 round(max(best_data)/10)*10],...
                        'xtick',1:10,'xticklabel',Info.blockKey(1:10))
                    titlestr = sprintf('%s | %s | ch %i, clu %i\n%i dBSPL | %i-%i Hz\nMTF RCorr, best filter sigma %i ms',...
                        subject,session,channel,clu,spl,SoundData(5,bkStart_ms(1)),lpn,best_sw);
                    title(titlestr)
                    xlabel('Block')
                    ylabel('Percent correctly assigned')
                    
                    
                    % Plot overall percent correct for each sigma value
                    hf2=figure;
                    bar(mean(percentCorrectMat(:,2,:),3),'k')
                    set(gca,'xticklabel',round(mean(percentCorrectMat(:,1,:),3)),...
                        'ylim',[0 ceil(max(mean(percentCorrectMat(:,2,:),3))/10)*10])
                    xlabel('sigma for gaussian smoothing filter')
                    ylabel('Overall percentage of blocks correctly assigned')
                    title(sprintf('%s | %s | ch %i, clu %i\n%i dBSPL | %i-%i Hz',...
                            subject,session,channel,clu,spl,SoundData(5,bkStart_ms(1)),lpn))
                    
                        
                    % Save figures
                    savename = sprintf('%s_%s_ch%i_clu%i_%idB_LP%ihz_RCorrMTF_sw%i',...
                        subject,session,channel,clu,spl,lpn,best_sw);
                    set(hf1,'PaperOrientation','landscape');
                    print(hf1,'-dpdf',fullfile(savedir,savename),'-bestfit')

                    savename = sprintf('%s_%s_ch%i_clu%i_%idB_LP%ihz_RCorr_overallPC',...
                        subject,session,channel,clu,spl,lpn);
                    set(hf2,'PaperOrientation','landscape');
                    print(hf2,'-dpdf',fullfile(savedir,savename),'-bestfit')

                    
                end %amd
            end %lpn
        end %spl
        
    end %clu
end %channel

end % sessions
