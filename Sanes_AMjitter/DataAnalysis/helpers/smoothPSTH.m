function [FRsmooth,peakFR,peakTime,minFR] = smoothPSTH(resp_jit)
% Input is mean number of spikes in 1 ms bins

FilterType =  'boxcar'; 'lowess'; 'gauss';
filtbinsize = 11; %odd number

FRraw = resp_jit.*1000;

%     figure(100); clf; hold on
%     plot(FRraw,'k')

if strcmp(FilterType,'boxcar')
    
    coeffFilt   = ones(1, filtbinsize)/filtbinsize;
    
    FRsmooth = filter(coeffFilt, 1, [FRraw zeros(1,20)]);
    
    fDelay    = (length(coeffFilt)-1)/2;
    adjustedT = [1:length(FRsmooth)]-fDelay;
    
    
    FRsmooth  = FRsmooth(adjustedT>0 & adjustedT<251);
    adjustedT = adjustedT(adjustedT>0 & adjustedT<251);
        
    [peakFR,peakTime] = max(FRsmooth);
    peakTime          = adjustedT(peakTime);
    minFR             = min(FRsmooth);
    
%         plot(adjustedT,FRsmooth,'r','LineWidth',2)
%         plot(peakTime,peakFR,'p','Color',[0.4 0 0],'MarkerSize',15,'LineWidth',2)
    
    
elseif strcmp(FilterType,'gauss')
    
    h = [1/2 1/2];
    binomialCoeff = conv(h,h);
    for n = 1:19
        binomialCoeff = conv(binomialCoeff,h);
    end
    
    FRsmooth = filter(binomialCoeff, 1, [FRraw zeros(1,20)]);
    
    fDelay    = (length(binomialCoeff)-1)/2;
    adjustedT = [1:length(FRsmooth)]-fDelay;
    
    FRsmooth  = FRsmooth(adjustedT>0 & adjustedT<251);
    adjustedT = adjustedT(adjustedT>0 & adjustedT<251);    
    
    [peakFR,peakTime] = max(FRsmooth);
    peakTime = adjustedT(peakTime);
    minFR             = min(FRsmooth);

    
%         plot(adjustedT,FRsmooth,'Color',[0 0 0.7],'LineWidth',2)
%         plot(peakTime,peakFR,'p','Color',[0 0 0.4],'MarkerSize',15,'LineWidth',2)
%         xlim([0 250])
    
        
elseif strcmp(FilterType,'lowess')
    
    FRsmooth = smooth(FRraw',filtbinsize,'loess')';
    
    [peakFR,peakTime] = max(FRsmooth);
     minFR            = min(FRsmooth);    
    
%         plot(1:250,FRsmooth,'Color',[0 0 0],'LineWidth',2)
%         plot(peakTime,peakFR,'p','Color',[0 0 0.4],'MarkerSize',15,'LineWidth',2)
%         xlim([0 250])
        
end


if size(FRsmooth,2) ~= 250
    keyboard
end
    
end




 