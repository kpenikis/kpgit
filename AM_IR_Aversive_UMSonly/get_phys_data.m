function filtPhys = get_phys_data(BLOCK,subject,tetN)
% called by: checkSortingQuality.m (which is still in progress)

global fn


%% Load this block epData datafile

block_str = sprintf('Block-%i.mat',BLOCK);
datafile = fullfile(fn.raw,subject,block_str);      %change back to fn.raw
if ~exist(datafile,'file')
    warning('raw epData file not found.')
    keyboard
end
fprintf(' loading data file %s...',datafile)
clear epData;
load(datafile,'-mat'); %loads data struct: epData
if isempty(epData)
    error('data file did not load correctly!')
else
    disp('done.')
end



keyboard

%% Only keep phys data from one tetrode

switch tetN
    case 1
        chs = 12:15;
    case 2
        chs = 8:11;
    case 3
        chs = 4:7;
    case 4
        chs = 1:3;
end

data_raw = epData.streams.Wave.data(chs,:);
fs = epData.streams.Wave.fs;
% clear epData;


%% Filter data (300 - 20k)


Wp = [ 300  6000] * 2 / fs;        %cutoff fqs (passband)
Ws = [ 225  8000] * 2 / fs;        %cutoff fqs (stopband)
[N,Wn] = buttord( Wp, Ws, 3, 20);  %create filter parameters
[B,A] = butter(N,Wn);              %build filter
fprintf('   filtering data...')
for ch = 1:size(data_raw,1)
    filtPhys(ch,:) = filtfilt( B, A, data_raw(ch,:) );
end


keyboard



end
