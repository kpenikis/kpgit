function [Stream_FRsmooth,Stream_zscore,Stream_Spks,ymaxval] = convertSpiketimesToFR(spiketimes,StreamNSamples,msStart,msSilEnd,SmoothType,SmoothWin,REFERENCE)
% [Stream_FRsmooth,Stream_zscore,Stream_Spks,ymaxval] = convertSpiketimesToFR(spiketimes,StreamNSamples,msStart,msSilEnd,bs_rect,bs_gaus,REFERENCE)
%
% For use with the long stream stimulus. Bins and smooths spiketimes to
% extract continuous firing rate throughout the session. Also calculates
% the z-score of FR, with reference to either the activity during the
% silent period at the beginning of each session, or relative to activity
% throughout the whole session. 
% Lastly, outputs the ymaxval, for plotting all rasters on the same axes.
%
% KP, 2018-03; updated 2019-12


if ~ischar(SmoothType)
    error('Update smoothing input variables')
end

Stream_Spks = zeros(1,StreamNSamples);
Stream_Spks(spiketimes) = 1;


% Get smoothed FR from spiketimes for entire stream
switch SmoothType
    
    case 'gauss'
        %gaussian window
        window=gausswin(SmoothWin);
        window=window-min(window);
        window=window/sum(window);
        
        Stream_FRsmooth = conv(Stream_Spks,window,'same').*1e3;
        
    case 'rect'
        %rectangular window
        window=rectwin(SmoothWin);
        window=window/sum(window);
        
        Stream_FRsmooth = conv(Stream_Spks,window,'same').*1e3;
        
    case 'exp'
        %exponential decay
        winlen = 1000;
        lambda = 1/SmoothWin;
        window = exp(-lambda*(1:winlen));
        
        Stream_FRsmooth = conv(Stream_Spks,window).*1e3;
        Stream_FRsmooth = Stream_FRsmooth(1:length(Stream_Spks));
end


% Convert FR to z-score

switch REFERENCE
    
    case 'stream'
        Stream_zscore = zscore(Stream_FRsmooth(msStart:end));
        Stream_zscore = [nan(1,msStart-1) Stream_zscore];
        
    case 'silence'
        
        % Get samples (ms) of silence
        ref_samps  = msStart:msSilEnd;
        
        % Calculate FR distribution for silent period in beginning
        meanFR = mean(Stream_FRsmooth(ref_samps));
        stdFR = std(Stream_FRsmooth(ref_samps));
        
%         if stdFR==0
%             stdFR = 0.1;
%             if meanFR==0
%                 meanFR = 0.1;
%             end
%         end
        
        Stream_zscore = (Stream_FRsmooth(msStart:end) - meanFR) / stdFR;
        Stream_zscore = [nan(1,msStart-1) Stream_zscore];
        
        %                 figure; plot(Stream_zscore,'k')
        %                 hold on
        %                 plot([TrialData.onset(2) TrialData.onset(2)],[-2 10],'g')
end


% Truncate Streams
Stream_Spks      = Stream_Spks(1:StreamNSamples);
Stream_FRsmooth  = Stream_FRsmooth(1:StreamNSamples);
Stream_zscore    = Stream_zscore(1:StreamNSamples);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set ymax based on overall FR of unit

yvals = [30 50 75 100 150 200 250 300 400 500 600 700 800 1000 1500];
yin = find( ( max(yvals, median(Stream_FRsmooth(msStart:end))+3*std(Stream_FRsmooth(msStart:end)) ) - yvals )==0 );
ymaxval = yvals(yin(1));



end
