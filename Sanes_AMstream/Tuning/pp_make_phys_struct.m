function ph = pp_make_phys_struct( epData, t1, t2, epoc_trs )
%
%  pp_make_phys_struct( epData, t1, t2 )
%    Sub-routine called by pp_prepare_format.
%    
%    Creates matrix of filtered ephys data. First, filters the entire data
%    stream from epData. Then extracts a window of data around each trial
%    onset.
%    
%
%  KP, 2016-04; last updated 2016-04
% 


% Get basic timing info

fs = epData.streams.Wave.fs;
wave_raw = double(epData.streams.Wave.data);


efns = fieldnames(epData.epocs);
onset_str = ['epData.epocs.' efns{1} '.onset'];
try
    onsets = round(fs*eval(onset_str)); %samples
    onsets = onsets(epoc_trs);
catch
    warning('Could not get trial start times.')
    keyboard
end



% Filter data

Wp = [ 300  6000] * 2 / fs;        %cutoff fqs (passband)
Ws = [ 225  8000] * 2 / fs;        %cutoff fqs (stopband)
[N,Wn] = buttord( Wp, Ws, 3, 20);  %create filter parameters
[B,A] = butter(N,Wn);              %build filter
fprintf('   filtering data...')
for ch = 1:size(wave_raw,1)
    if ch==8;    continue;      end
    wave_filt(ch,:) = filtfilt( B, A, wave_raw(ch,:) );
end



%%% For UMS, data structure must be either:
%%%    matrix format [trials x samples x channels]
%%%    or cell array {trials}[samples x channels]

ph = nan( numel(onsets), 1+t2-t1, size(wave_filt,1) );

for it = 1:numel(onsets)  %iterate through trials
    
    t0 = onsets(it);
    
    for ic = 1:size(wave_filt,1)  %iterate through channels

        try
            if (t0+t2) > size(wave_filt,2)  %if the recording ends before the extraction window
%                 keyboard
                ph(it, 1 : (1+t2-t1-((t0+t2)-size(wave_filt,2))) ,ic) = wave_filt(ic, (t0+t1):end );
                
            else
                % Pull window of data around trial onset
                ph(it,:,ic) = wave_filt(ic, (t0+t1):(t0+t2) );
                
            end
            
        catch
            warning('inconsistent window sizes')
            keyboard
            
        end
        
    end
end


if any(isnan(ph(it,:,ic)))
    warning('recording cut off early. mark last trial for removal')
%     keyboard
%     ph = ph(1:(it-1),:,:);  % must you remove the last row of Stim struct too then??
end


end
