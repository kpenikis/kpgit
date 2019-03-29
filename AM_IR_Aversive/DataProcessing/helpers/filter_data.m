function dataFilt = filter_data( data, fs )
%  DOESNT WORK BECAUSE DATA MATRIX TOO BIG. WRITING FILTERED DAT FILE
%  VIA KILOSORT INSTEAD.
% 



% Filter data
dataFilt = nan(size(data));

Wp = [ 300  6000] * 2 / fs;        %cutoff fqs (passband)
Ws = [ 225  8000] * 2 / fs;        %cutoff fqs (stopband)
[N,Wn] = buttord( Wp, Ws, 3, 20);  %create filter parameters
[B,A] = butter(N,Wn);              %build filter

for ch = 1:size(data,1)
    dataFilt(ch,:) = filtfilt( B, A, data(ch,:) );
end



end
