

whichStim  = 'AC';
whichClass = 'SumAV';

% Set N trials
allMT = [12 16 17 19 22 32];

for minTrs = allMT
close all

% Run MC_eachCell
% MC_eachCell(minTrs,whichStim,whichClass)
% close all

% Run MC_subpop
MC_subpop(minTrs,whichStim,whichClass)

end
