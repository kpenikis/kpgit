function rcorr_class_mph(MPH,subject,session,channel,clu,spl,lpn,AMrate)



% Determine which stimuli to compare
% (time in block 0 or 1000 and context)

T = cell(2,length(blocksN));
context = fieldnames(MPH);
data = struct;
ii = 0;
for icontext = context'
    
    % Check for early first
    if sum( MPH.(icontext{:}).pdc.pdtime(:,2) == 0 )>12 &&...
            all((MPH.(icontext{:}).pdtime(MPH.(icontext{:}).pdtime(:,1) == 1,2) - 0)<100)
        ii=ii+1;
        data(ii).raster = 
    end
    
    
end


keyboard


for icontext = context'
    for ii = 1:numel(unique(MPH(ir).(icontext{:}).pdtime(:,1)))
        
        these_rows = find(MPH(ir).(icontext{:}).pdtime(:,1)==ii)';
        
        [~,sprow] = min(abs( pdtimes - round(mode(MPH(ir).(icontext{:}).pdtime(these_rows,2))) ));
        sp_idx = 2*sprow + 2;
        
        %~~~~~~~~~~
        % IR first
        %~~~~~~~~~~
        % Get raster data
        raster_x = []; raster_y = [];
        for it = 1:numel(these_rows)
            raster_x = [raster_x find( MPH(ir).(icontext{:}).raster(these_rows(it),:) )];
            raster_y = [raster_y repmat(nrows+it,1,sum( MPH(ir).(icontext{:}).raster(these_rows(it),:) ))];
        end
        
        % Plot rasters (left early, right middle)
        ax_r(2)=subplot(9,2,2);
        hold on
        hL=plot(raster_x, raster_y,'+','Color',contextcolors(1+find(strcmp(icontext,context)),:),...
            'MarkerSize',2,'LineWidth',0.5);
        % If first column, set axis labels
        if ii==1
            xlabel('time (ms)')
            ylabel([num2str(min(blocks_N)) ' trs each'])
        end
        %                                 % Set marker transparency
        %                                 if ~isempty(hL)
        %                                     hMarkers=[];
        %                                     while isempty(hMarkers)
        %                                         hMarkers = hL.MarkerHandle;
        %                                         pause(0.1)
        %                                     end
        %                                     hMarkers.EdgeColorData(4) = uint8(255 * 0.4 );
        %                                 end
        
        % Plot folded MPH
        ax_m(sp_idx-2) = subplot(9,2,sp_idx);
        sphist = reshape([hist(raster_x,linspace(0,1000/AMrates(ir),52)); hist(raster_x,linspace(0,1000/AMrates(ir),52))],1,52*2);
        patch(reshape([0.5:52.5; 0.5:52.5],1,53*2), [0 sphist 0], contextcolors(1+find(strcmp(icontext,context)),:) )
        %                                 bar(hist(raster_x,linspace(0,1000/AMrates(ir),52)),1,...
        %                                     'FaceColor',contextcolors(1+find(strcmp(icontext,context)),:),...
        %                                     'EdgeColor',contextcolors(1+find(strcmp(icontext,context)),:))
        box off
        %                                 ax_m(sp_idx-2).YLim(1) = 0;
        title([icontext{:} ', ' num2str(round(mode(MPH(ir).(icontext{:}).pdtime(these_rows,2)))) 'ms'])
        
        % Keep track of largest bin spk count to
        % set ymax
        if max(hist(raster_x,linspace(0,1000/AMrates(ir),52))) > max_count
            max_count = max(hist(raster_x,linspace(0,1000/AMrates(ir),52)));
        end
        
        
        %~~~~~~~~~~
        % PDC now
        %~~~~~~~~~~
        
        these_rows = find(MPH(ir).(icontext{:}).pdc.pdtime(:,1)==ii)';
        
        % Get raster data
        raster_x = []; raster_y = [];
        for it = 1:numel(these_rows)
            raster_x = [raster_x find( MPH(ir).(icontext{:}).pdc.raster(these_rows(it),:) )];
            raster_y = [raster_y repmat(nrows+it,1,sum( MPH(ir).(icontext{:}).pdc.raster(these_rows(it),:) ))];
        end
        
        % Plot rasters (left early, right middle)
        ax_r(1)=subplot(9,2,1);
        hold on
        hL=plot(raster_x, raster_y,'+','Color',contextcolors(1,:),...
            'MarkerSize',2,'LineWidth',0.5);
        % If first column, set axis labels
        if ii==1
            xlabel('time (ms)')
            ylabel([num2str(min(blocks_N)) ' trs each'])
        end
        %                                 % Set marker transparency
        %                                 if ~isempty(hL)
        %                                     hMarkers=[];
        %                                     while isempty(hMarkers)
        %                                         hMarkers = hL.MarkerHandle;
        %                                         pause(0.1)
        %                                     end
        %                                     hMarkers.EdgeColorData(4) = uint8(255 * 0.4 );
        %                                 end
        
        % Plot folded MPH
        ax_m(sp_idx-1-2) = subplot(9,2,sp_idx-1);
        sphist = reshape([hist(raster_x,linspace(0,1000/AMrates(ir),52)); hist(raster_x,linspace(0,1000/AMrates(ir),52))],1,52*2);
        patch(reshape([0.5:52.5; 0.5:52.5],1,53*2), [0 sphist 0], contextcolors(1,:) )
        %                                 bar(hist(raster_x,linspace(0,1000/AMrates(ir),52)),1,...
        %                                     'FaceColor',contextcolors(1,:),...
        %                                     'EdgeColor','k')
        box off
        %                                 ax_m(sp_idx-2).YLim(1) = 0;
        title(['periodic, ' num2str(round(mode(MPH(ir).(icontext{:}).pdc.pdtime(these_rows,2)))) 'ms'])
        
        % Keep track of largest bin spk count to
        % set ymax
        if max(hist(raster_x,linspace(0,1000/AMrates(ir),52))) > max_count
            max_count = max(hist(raster_x,linspace(0,1000/AMrates(ir),52)));
        end
        
        
        
        % Add to count of total number of rows
        nrows = nrows+min(blocks_N);
        
    end %ii
    
    
end %context


%~~~~~~~~
sws = 16;
nTemps = 10;
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



end