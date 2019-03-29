function reset_ch_sort(subject,session,channel)


% Load data 
fn = set_paths_directories;
fprintf('loading data...\n')
if ~exist('Spikes','var')
    filename = sprintf('%s_sess-%s_Spikes',subject,session);
    load(fullfile(fn.processed,subject,filename));
end
filename = sprintf('%s_sess-%s_Info',subject,session);
load(fullfile(fn.processed,subject,filename));


% Reset sort code
Spikes.man_sort(channel) = 0; 


% Save Spikes structure
fprintf('\nsaving data...\n')
savename = sprintf('%s_sess-%s_Spikes',subject,session);
save(fullfile(fn.processed,subject,savename),'Spikes','-v7.3');

end
