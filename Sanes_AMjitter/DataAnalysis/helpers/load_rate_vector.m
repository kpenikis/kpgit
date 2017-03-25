
function rV = load_rate_vector(blocks,subject)

stimdir = fullfile('/Users/kpenikis/Documents/SanesLab/Data/AMJitter/RawData',subject,sprintf('Block-%i_Stim',blocks(1)));
rV = load(fullfile(stimdir,stim(1).stimfn));


