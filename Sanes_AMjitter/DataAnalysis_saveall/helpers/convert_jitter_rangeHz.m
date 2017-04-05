function rangeHz = convert_jitter_rangeHz(jitter_ID,center_rate_Hz)
% rangeHz = convert_jitter_rangeHz(jitter_ID,center_rate_Hz)
%   takes integers that were used to ID the jitter vectors, and converts
%   them to the corresponding range of AM rates in each vector.

rangeHz = 2.^((log2(center_rate_Hz))+jitter_ID./100) - 2.^((log2(center_rate_Hz))-jitter_ID./100);

end