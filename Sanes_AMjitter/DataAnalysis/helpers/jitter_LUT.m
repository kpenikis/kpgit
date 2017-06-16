function rateVector = jitter_LUT(rate,jit_str)
% Get rate vector for specified stimulus
%  jitter_LUT(CenterRate,jit_str)
%   KP 04-27-2017

fn = set_paths_directories;

% Find stimulus buffer file
stfn = ['AM_' num2str(rate) 'Hz_' jit_str '.mat'];
if ~exist(fullfile(fn.stim,stfn),'file')
    warning('cannot find stimulus file in directory')
else
    rV = load(fullfile(fn.stim,stfn));
end

rateVector = rV.buffer(2:end);

end

