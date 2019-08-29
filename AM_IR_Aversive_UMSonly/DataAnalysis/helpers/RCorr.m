function [best_data,best_PC,best_sw] = RCorr(StreamSpikes,AllStimStarttimes,minDur,nIterations,subject,session,channel,clu)
%
%  ap_rcorr(subject)
%
%
%  KP, 2018-03
%

global fn 

% Set some plotting params
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');  % [ left bottom width height ]
figsize = [1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)/2];


%~~~~~~~~
sws = 16; %[1 2 4 8 16 32 64 128 256];
nTemps = nIterations;
mspad  = 500;
%~~~~~~~~

percentCorrectMat = nan(numel(sws),2,nTemps);
best_PC   = 0;
best_sw   = 0;
best_data = [];

nstim = numel(AllStimStarttimes);

for isw = 1:numel(sws)
    
    % Create gaussian smoothing window for convolving spiketrains
    sw = sws(isw);
    [~,GW,~] = rc_makeSmoothGauss(sw,minDur);
    
    % Determine the half length (in samples) of the gaussian window. Use later
    % for shifting convolved signal to align with correct time points.
    lGausHalf = (length(GW)-1)/2;
    
    % Create empty matrix to track the template trials used
    ntUsed = zeros(1,nstim);
    
    percentAssignment = nan(nstim,nstim,nTemps);
    
    h = waitbar(0,sprintf('cycling through templates for GW=%i...',sw));
    
    for iTemp = 1:nTemps
        
        waitbar(iTemp/nTemps,h)
        
        [T,ntUsed] = rc_getTemplates(StreamSpikes,GW,AllStimStarttimes,minDur,ntUsed,mspad);
        
        
        % Preallocations
        finalMat = [];
        count = 0;
        
        % Go through each stimulus
        for is = 1:nstim
            
            Starttimes = AllStimStarttimes{1,is};
            if isempty(Starttimes)
                continue
            end
            
            % Go through each trial
            for it = 1:numel(Starttimes)
                
                % Skip the current template trial
                if it==ntUsed(end,is)
                    continue
                end
                
                % Get spiketimes
                bkStart_ms = Starttimes(it);
                bkStop_ms = bkStart_ms+minDur-1;
                
                % Get spiketimes during the selected trial, adjusted to start time 0
                sp01=[];
                sp01 = StreamSpikes((bkStart_ms-mspad):(bkStop_ms+mspad));
                
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
                
                [blockAssignment,maxR] = rc_calcR(T,S);
                count = count+1;
                finalMat(count,:) = [is,it,blockAssignment,maxR];
                
                
            end %it
            
            
        end %is
        
        
        % Calculate overall percent correct
        correct = find(finalMat(:,1) == finalMat(:,3));
        percentCorrect = (numel(correct)/length(finalMat))*100;
        percentCorrectMat(isw,:,iTemp) = [sw,percentCorrect];
        
        
        % Collect data to visualize full results
        for thisstim = 1:max(finalMat(:,1))
            percentAssignment(thisstim,:,iTemp) = hist(finalMat(finalMat(:,1) == thisstim,3),1:nstim)./sum(finalMat(:,1) == thisstim);
        end
        
        
    end %iTemp
    
    
    % If this sw has a higher overall PC than current,
    % save the data to plot later in MTF
    if mean(percentCorrectMat(isw,2,:),3)>best_PC
        
        best_sw   = sw;
        best_PC   = mean(percentCorrectMat(isw,2,:),3);
        best_data = diag(mean(percentAssignment,3)).*100;
        
% %         % Plot full results, averaged across the random
% %         % template trials chosen
% %         hf=figure; hold off
% %         imagesc(mean(percentAssignment,3)')
% %         axis square
% %         set(gca,'xtick',1:nstim,'ytick',1:nstim)
% %         titlestr = sprintf('%s | %s | ch %i, clu %i\nRCorr results, with filter sigma %i ms',...
% %             subject,session,channel,clu,sw);
% %         title(titlestr)
% %         xlabel('True stimulus')
% %         ylabel('Assigned stimulus')
% %         colormap('bone')
% %         colorbar
% %         caxis([0 0.5])
% %         
% %         % Save figure
% %         savedir = fullfile(fn.processed,'RCorr');
% %         if ~exist(savedir,'dir')
% %             mkdir(savedir)
% %         end
% %         savename = sprintf('%s_%s_ch%i_clu%i_RCorrMat_sw%i',...
% %             subject,session,channel,clu,sw);
% %         print_eps_kp(hf,fullfile(savedir,savename))
        
    end
    
    close(h)
    
end %isw



% [maxPC,imaxPC] = max(best_data);

% % And plot the MTF/tuning function of the best sw
% hf1=figure; hold on
% set(hf1,'Position',figsize)
% plot(1,best_data(9),'k.--','MarkerSize',30,'LineWidth',1)
% plot(2,best_data(1),'k.--','MarkerSize',30,'LineWidth',1)
% plot(3:7,best_data(2:6),'k.--','MarkerSize',30,'LineWidth',1)
% plot(8:9,best_data(7:8),'k.--','MarkerSize',30,'LineWidth',1)
% set(gca,'xlim',[0 10],'ylim',[0 ceil(max(best_data)/10)*10],...
%     'xtick',1:9,'xticklabel',{'Sil' 'Warn' '2 Hz' '4 Hz' '8 Hz' '16 Hz' '32 Hz' 'AC' 'DB'})
% titlestr = sprintf('%s | %s | ch %i, clu %i\nMTF RCorr, best filter sigma %i ms',...
%     subject,session,channel,clu,best_sw);
% title(titlestr)
% xlabel('Stimulus')
% ylabel('Percent correctly assigned')
% 
% 
% % Plot overall percent correct for each sigma value
% hf2=figure;
% bar(mean(percentCorrectMat(:,2,:),3),'k')
% set(gca,'xticklabel',round(mean(percentCorrectMat(:,1,:),3)),...
%     'ylim',[0 ceil(max(mean(percentCorrectMat(:,2,:),3))/10)*10])
% xlabel('sigma for gaussian smoothing filter')
% ylabel('Overall percentage of blocks correctly assigned')
% title(sprintf('%s | %s | ch %i, clu %i',...
%     subject,session,channel,clu))
% 
% 
% % Save figures
% savename = sprintf('%s_%s_ch%i_clu%i_RCorrMTF_sw%i',...
%     subject,session,channel,clu,best_sw);
% set(hf1,'PaperOrientation','landscape');
% print(hf1,'-dpdf',fullfile(savedir,savename),'-bestfit')
% 
% savename = sprintf('%s_%s_ch%i_clu%i_RCorr_overallPC',...
%     subject,session,channel,clu);
% set(hf2,'PaperOrientation','landscape');
% print(hf2,'-dpdf',fullfile(savedir,savename),'-bestfit')


end %function
