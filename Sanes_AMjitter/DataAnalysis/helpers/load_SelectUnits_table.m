function SelectUnits = load_SelectUnits_table(SU,var1,threshold)

fn = set_paths_directories;

if SU
    savedir = fullfile(fn.standardPd,'onlySU','comparison',var1);
else
    savedir = fullfile(fn.standardPd,'comparison',var1);
end

tablesavename = sprintf('SelectUnits_%s_%i',var1,threshold);

SelectUnits = readtable(fullfile(savedir,tablesavename));


end






