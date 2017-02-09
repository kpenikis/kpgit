function cludata = ap_nmData_depth(cludata)

% Plot 1) classifier PC fit and d' transformed, 1b) classifier d' formula fit, 2) FR d' formula fit.

% cludata.nmdata.classifier
% cludata.nmdata.formula

%   ~~ jitters, timebins ~~
% cludata.nmdata.classifier.output.PCdata{j}  =  [  stimval  #yesAM  #trials  ]
% cludata.nmdata.classifier.output.dprime{j}  =  [ stimval X d' ]
% cludata.nmdata.classifier.output.stim{j,1}  =  'jitter'  
% cludata.nmdata.classifier.output.stim{j,2}  =  [depths]

% cludata.nmdata.formula.FRmat{j}
% cludata.nmdata.formula.dpmat{j}


%% CLASSIFIER DATA
metric   = {'FR' 'SpV' 'SpV' 'SpV' 'SpV' 'SpV'};
binsizes = [ 0    250   100   50    25    10];

for ib = 1:numel(binsizes)
    
    thisbin = binsizes(ib);
    thismet = metric{ib};
    
    % Run classifier for all stimuli of this session
    fprintf('   Running classifier for %s binsize %i...\n',thismet,thisbin)
    
    output = get_classifier_data( cludata.raster, thismet, thisbin, 0);
    
    cludata.nmdata.classifier(ib).params.metric  = thismet;
    cludata.nmdata.classifier(ib).params.binsize = thisbin;
    cludata.nmdata.classifier(ib).output = output;
    
end

fprintf('#Finished.\n\n')


%% FORMULA DATA







end