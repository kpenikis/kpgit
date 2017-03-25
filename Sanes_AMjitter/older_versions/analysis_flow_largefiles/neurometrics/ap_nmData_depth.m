function [nmdata,cludata] = ap_nmData_depth(cludata,is,ip)

% Plot 1) classifier PC fit and d' transformed, 1b) classifier d' formula fit, 2) FR d' formula fit.


%% CLASSIFIER DATA
nmdata = struct();

metric   = {'FR' 'SpV' 'SpV' 'SpV' 'SpV' 'SpV'};
binsizes = [ 0    250   100   50    25    10];

for ib = 1:numel(binsizes)
    
    thisbin = binsizes(ib);
    thismet = metric{ib};
    
    fprintf('   Running classifier for %s binsize %i...\n',thismet,thisbin)
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [this_output,cludata] = format_classifier_data( cludata, is, ip, thismet, thisbin );
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    nmdata.classifier(ib).params.metric  = thismet;
    nmdata.classifier(ib).params.binsize = thisbin;
    nmdata.classifier(ib).output = [this_output.output];
    
end



%% FORMULA DATA

fprintf('   Getting FR formula data...\n')

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[this_output,cludata] = format_formula_data( cludata, is, ip );
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nmdata.formulaFR.output = [this_output.output];


fprintf('#Finished.\n\n')
















end