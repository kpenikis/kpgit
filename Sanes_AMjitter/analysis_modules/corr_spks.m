function [Rs,ps,shift,trR,Waves] = corr_spks(stim,subject)

% Get stimulus waveforms (fs: 1000 Hz)
Waves = ap_stimplotting(subject,stim);

% Set up empty vectors for data
shift = -200:5:200;
Rs = nan( numel(stim), numel(shift) );
ps = nan( numel(stim), numel(shift) );

% Go through each stimulus
for is = 1:numel(stim)
    
    data = stim(is);
    wave = Waves(is);
    

    if data.AMdepth==0
        wave.y = wave.y + (0.0002.*rand(size(wave.y))-0.0001);
    end
    
    for sh = shift
    
    r = nan(max(data.y),1);
    p = nan(max(data.y),1);
    for it = 1:max(data.y)
        x = data.x(data.y==it) - sh; 
        x = x(x>0 & x<=wave.x_ms(2));
        
        spBinary(:,it) = zeros(wave.x_ms(2),1);
        spBinary(x,it) = 1;
        
        [r(it),p(it)] = corr(spBinary(:,it),wave.y');
        
    end
    
    Rs(is,sh==shift) = nanmean(r);
    ps(is,sh==shift) = nanmean(p);
    
    if sh==0
        fullR = corr(spBinary);
        trR(is,1) = mean(fullR(fullR~=1));
    end
    
    clear r p spBinary
    end
    
end


end


