
MergedUnits = struct;

iu = 0;

%
newSessLabel = 'SX';
AcceptMerge = 0;
AcceptMerge = checkUnits( 'WWWf_253400', {'SB' 'SC'}, 3, [17 45], newSessLabel);
if AcceptMerge
    iu = iu+1;
    MergedUnits(iu).new = newSessLabel;
    MergedUnits(iu).un1 = 'WWWf_SB_3_17';
    MergedUnits(iu).un2 = 'WWWf_SC_3_45';
end

%
newSessLabel = 'YX';
AcceptMerge = 0;
AcceptMerge = checkUnits( 'WWWf_253400', {'YA' 'YB'}, 11, [11 5], newSessLabel);
if AcceptMerge
    iu = iu+1;
    MergedUnits(iu).new = newSessLabel;
    MergedUnits(iu).un1 = 'WWWf_YA_11_11';
    MergedUnits(iu).un2 = 'WWWf_YB_11_5';
end

%
newSessLabel = 'BX';
AcceptMerge = 0;
AcceptMerge = checkUnits( 'AAB_265055', {'BA' 'BB'}, 4, [32 23], newSessLabel);
if AcceptMerge
    iu = iu+1;
    MergedUnits(iu).new = newSessLabel;
    MergedUnits(iu).un1 = 'AAB_BA_4_32';
    MergedUnits(iu).un2 = 'AAB_BB_4_23';
end

%
newSessLabel = 'BX';
AcceptMerge = 0;
AcceptMerge = checkUnits( 'AAB_265055', {'BA' 'BB'}, 11, [16 1], newSessLabel);
if AcceptMerge
    iu = iu+1;
    MergedUnits(iu).new = newSessLabel;
    MergedUnits(iu).un1 = 'AAB_BA_11_16';
    MergedUnits(iu).un2 = 'AAB_BB_11_1';
end

fn = set_paths_directories;
save(fullfile(fn.processed,'MergedUnits'),'MergedUnits','-v7.3')

SplitUnits = MergedUnits;
save(fullfile(fn.processed,'SplitUnits'),'SplitUnits','-v7.3')



% subjects = { 'AAB_265055' 'AAB_265059' 'WWWf_253400' 'WWWlf_253395'  };
