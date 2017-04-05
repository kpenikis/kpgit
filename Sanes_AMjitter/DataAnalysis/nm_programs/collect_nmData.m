function nmdata = collect_nmData(rasterdata,indVar)
% called by getALL_nmData
% Plot 1) classifier PC fit and d' transformed, 1b) classifier d' formula fit, 2) FR d' formula fit.


%% CLASSIFIER DATA
nmdata = struct();

% metric   = {'FR' 'SpV' 'SpV' 'SpV' 'SpV' 'SpV'};
% binsizes = [ 0    250   100   50    25    10];

metric   = {'FR' 'SpV' 'SpV' 'SpV' 'Corr'};
binsizes = [ 0    100    25    10     0  ];

% metric   = {'Corr'};
% binsizes = [  0  ];

for ib = 1:numel(binsizes)
    
    thisbin = binsizes(ib);
    thismet = metric{ib};
    
    fprintf('   Running classifier for %s binsize %i...\n',thismet,thisbin)
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [this_output,stimC] = format_classifier_data( rasterdata, thismet, thisbin, indVar );
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nmdata.stim = stimC;
    fldnm = sprintf('classif%s%i',thismet,thisbin);
    nmdata.(fldnm) = [this_output.output];
    
end



%% FORMULA DATA

fprintf('   Getting FR formula data...\n')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[this_output,stimF] = format_formula_data( rasterdata, indVar );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nmdata.formulaFR = [this_output.output];


% Compare the stim data
cmprstim = cellfun(@strcmp,stimC,stimF,'UniformOutput',false);
if any([cmprstim{1}; cmprstim{2}] ~=1)
    warning('stim cells may not be identical')
    keyboard
end




end


