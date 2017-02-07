
function raster = get_raster(subject,session,channel,clu)

datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
dataname = sprintf('%s_sess-%s_Data',subject,session);
load(fullfile(datadir,subject,dataname))

iclu = find((Data.ch(channel).clu.label(:,1))==clu);
raster = Data.ch(channel).clu(iclu).raster;

end