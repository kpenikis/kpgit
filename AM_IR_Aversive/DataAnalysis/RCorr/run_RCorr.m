function PCMat = run_RCorr(StreamSpikes,AllStimStarttimes)
%
%  called by 
%
%
%  KP, 2018-03
%

global nIterations StimDur exclOnset stids

% Set some plotting params
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)

scrsz = get(0,'ScreenSize');  % [ left bottom width height ]
figsize = [1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)/2];


%~~~~~~~~
sws = 16; %[1 2 4 8 16 32 64 128 256]; %debug before running with >1 sws
mspad  = 500;
%~~~~~~~~

nstim = numel(AllStimStarttimes);

percentCorrectMat = nan(nstim,nstim,numel(sws),nIterations);
best_PC   = 0;
best_sw   = 0;
best_data = [];

for isw = 1:numel(sws)
    
    % Create gaussian smoothing window for convolving spiketrains
    sw = sws(isw);
    [~,GW,~] = rc_makeSmoothGauss(sw,StimDur);
    
    % Determine the half length (in samples) of the gaussian window. Use later
    % for shifting convolved signal to align with correct time points.
    lGausHalf = (length(GW)-1)/2;
    
    % Create empty matrix to track the template trials used
    ntUsed = zeros(1,nstim);
    
    percentAssignment = nan(nstim,nstim,nIterations);
    
    h = waitbar(0,sprintf('cycling through templates for GW=%i...',sw));
    
    for iIt = 1:nIterations
        
        waitbar(iIt/nIterations,h)
        
        [T,ntUsed] = rc_getTemplates(StreamSpikes,GW,AllStimStarttimes,StimDur,ntUsed,mspad);
        
        stids = find(~cellfun(@isempty,T(1,:)));
        
        % Go through each stimulus
        for isTrue = stids
            
            Starttimes = AllStimStarttimes{1,isTrue};
            if isempty(Starttimes)
                continue
            end
            
            % Preallocations
            finalMat = [];
            count = 0;
            
            % Go through each trial
            for it = 1:numel(Starttimes)
                
                % Skip the current template trial
                if it==ntUsed(end,isTrue)
                    continue
                end
                
                % Get spiketimes
                bkStart_ms = Starttimes(it);
                bkStop_ms = bkStart_ms+StimDur-1;
                if exclOnset
                    bkStart_ms = bkStart_ms+150;
                end
                
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
                S = sp_conv(mspad+(1:StimDur));
                
                % Check by plotting again
                %     figure; hold on
                %     h = plot(sp_conv,'b-');
                %     y = zeros(length(sp),1);
                %     h2 = plot(sp-mspad,y,'b.');
                %     xlim([1 StimDur])
                
                
                % Calculate RCORR values and make block assignment
                
                [blockAssignment,maxR] = rc_calcR(T,S);
                count = count+1;
                finalMat(count,:) = [it blockAssignment maxR];
                
                
                
            end %it
            
            %percentCorrectMat = nan(nstim,nstim,numel(sws),nIterations);
            for isAss =  1:nstim
                if ismember(isAss,stids)
                    percentCorrectMat(isTrue,isAss,isw,iIt) = sum(finalMat(:,2)==isAss)/size(finalMat,1);
                end
            end
            
        end %isTrue
        
    end %iIt
    
end %isw
close(h)

PCMat = mean(percentCorrectMat(:,:,1,:),4);


% figure; hold off
% imagesc(PCMat)
% caxis(gca,[0 1])
% colormap(gca,cmocean('algae'))
% colorbar
% xlabel('True Stim')
% ylabel('Assigned')
% axis square


end %function
