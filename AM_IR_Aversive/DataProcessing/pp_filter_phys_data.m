function ph = pp_filter_phys_data( epData )
%
%
%  KP, 2017-05; last updated 2016-04
% 


% Get basic timing info

fs = epData.streams.Wave.fs;
wave_raw = double(epData.streams.Wave.data);


% Filter data

Wp = [ 300  6000] * 2 / fs;        %cutoff fqs (passband)
Ws = [ 225  8000] * 2 / fs;        %cutoff fqs (stopband)
[N,Wn] = buttord( Wp, Ws, 3, 20);  %create filter parameters
[B,A] = butter(N,Wn);              %build filter
fprintf('   filtering data...')
for ch = 1:size(wave_raw,1)
%     if ch==16;    continue;      end
    ph(ch,:) = filtfilt( B, A, wave_raw(ch,:) );
end



end
