function load_ZdataRdata(subject)

% subject = 'WWWf_244303';
fn = set_paths_directories(subject);

savedir = fullfile(fn.processed,subject,'^an_plots','Population');
Zdata = readtable(fullfile(savedir,sprintf('zFR_ObsPred_SU_%s',subject)));

savedir = fullfile(fn.processed,subject,'^an_plots','MPH','corr');
Rdata = readtable(fullfile(savedir,sprintf('MPH-corr_pdcIR_%s',subject)));




Zdata.zFR_normdiff = (Zdata.zFR_obs-Zdata.zFR_pred)./abs(Zdata.zFR_pred);
sortrows(Zdata,'zFR_normdiff')


