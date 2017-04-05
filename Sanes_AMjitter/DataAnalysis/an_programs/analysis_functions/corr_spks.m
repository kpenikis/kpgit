function [Rs,ps,shift,R_eachtr,trtrR,Waves] = corr_spks(raster,subject,Lag_for_R_eachtr)
% [Rs,ps,shift,R_eachtr,trtrR,Waves] = corr_spks(stim,subject,Lag_for_R_eachtr)
% 
%  Takes a raster format input. Also needs string of subject name. Optional
%  input of the lag at which to compute the trial-by-trial correlation
%  value (as called by get_classifier_data).
%  
%  Outputs:
%   Rs: mean correlation values averaged across trials, at each shift lag. 
%   ps: mean corr significance values avg across trials, at each shift lag.
%   shift: shift/lag values used.
%   R_eachtr: correlation values at ONE set lag, for each trial.
%   trtrR: trial-to-trial correlation, independent of waveform.
%   Waves: structure containing info about stimulus waveform.
%   
%  KP 2017.
%

if nargin<3
    R_eachtr = nan;
end

% Get stimulus waveforms (fs: 1000 Hz)
Waves = ap_stimplotting(subject,raster);

% Set up empty vectors for data
shift = -200:5:200;
Rs = nan( numel(raster), numel(shift) );
ps = nan( numel(raster), numel(shift) );

% Go through each stimulus
for istim = 1:numel(raster)
    
    data = raster(istim);
    wave = Waves(istim);
    

    if data.AMdepth==0
        wave.y = wave.y + (0.00002.*rand(size(wave.y))-0.00001);
    end
    
    for sh = shift
    
    r = nan(max(data.y),1);
    p = nan(max(data.y),1);
    if nargin>2 && sh==Lag_for_R_eachtr
        R_eachtr = nan(max(data.y),1);
    end
    for it = 1:max(data.y)
        x = data.x(data.y==it) - sh; 
        x = x(x>0 & x<=wave.x_ms(2));
        
        spBinary(:,it) = zeros(wave.x_ms(2),1);
        spBinary(x,it) = 1;
        
        [r_out,p_out] = corrcoef(spBinary(:,it),wave.y');
        if numel(r_out)==4
            r(it) = r_out(2,1);
            p(it) = p_out(2,1);
        else
            keyboard
        end
        if nargin>2 && sh==Lag_for_R_eachtr
            R_eachtr(it) = r_out(2,1);
        end
        
    end
    
    % Mean correlation across trials for this shift
    Rs(istim,sh==shift) = mean(r,'omitnan');
    ps(istim,sh==shift) = mean(p,'omitnan');
    
    if sh==0
        fullR = corrcoef(spBinary);
        trtrR(istim,1) = mean(fullR(fullR~=1),'omitnan');
    end
    
    clear r p spBinary
    end
    
end


end


