function Data = ap_nmPlotsAll_depth(subject,session)

%~~~~~~~~~~~~~~~
FORCECLASS = 1;
%~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~
FORCEFIT   = 0;
%~~~~~~~~~~~~~~~

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers')
addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')


% Load Data file for session
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
filename = sprintf( '%s_sess-%s_Data',subject,session); 
load(fullfile(datadir,subject,filename))

%%  Calculate neurometric data 

% Step through each clu of each channel
for channel = 1:numel(Data.ch)
    for iclu = 1:numel(Data.ch(channel).clu)
        
        if isempty(Data.ch(channel).clu)
            continue
        end
        
        cludata = Data.ch(channel).clu(iclu);
        
        %  Calculate neurometric data if not yet run, or if rerun specified.
        if ~isfield(cludata,'nmdata') || FORCECLASS
            
            fprintf('Getting neurometric data for channel %i clu %i\n',channel,Data.ch(channel).clu(iclu).label(1))
            cludata = rmfield(cludata,'nmdata');
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            cludata = ap_nmData_depth(cludata);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Add to Data sturct and save file
            Data.ch(channel).clu(iclu).nmdata = cludata.nmdata;
            save(fullfile(datadir,subject,filename),'Data','-v7.3');
                        
        end
        
    end
end



%       ***   MUST ADD FORMULA DATA TOO   ***



%{

%%  Fit neurometric data 

[options, plotOptions] = setOptions;
options.dprimeThresh = 1;

% Step through each clu of each channel
for channel = 11%:numel(Data.ch)
    for iclu = 1:numel(Data.ch(channel).clu)
        
        if isempty(Data.ch(channel).clu)
            continue
        end
        
        for icp = 1:numel(Data.ch(channel).clu(iclu).nmdata.classifier)
            
            params = Data.ch(channel).clu(iclu).nmdata.classifier(icp).params;
            output = Data.ch(channel).clu(iclu).nmdata.classifier(icp).output;
            
            for is = 1:size(output.stim,1)
                
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Generate fit of Proportion Yes data
                [options,results,zFA] = fit_PYdata(output.PYdata{is},options);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
                %%   HERE BEGIN EDITS
                
                f1=figure; handles_f1=[];
                
                % Plot the percent correct values and fit
                figure(f1);
                s1 = subplot(2,2,1);
                plotPsych(results,plotOptions);
                handles_f1 = [handles_f1;s1];
                
                % Transform to dprime space and plot
                
                
                figure(f1)
                s2 = subplot(2,2,2);
                [x, fitted_yes, fitted_dprime, threshold, slope] =   ...
                    plotPsych_dprime( results, output.dprime{is}, ...
                                      options, plotOptions, zFA);
                hold on;
                
                %Save plot handles and link axes
                handles_f1 = [handles_f1;s2]; %#ok<*AGROW>
                linkaxes(handles_f1,'x');
                
                %Save figure
                suptitle('Classifier disriminability')
                fname = [filename(1:end-4),'_psychometric_fits'];
                set(f1,'PaperPositionMode','auto');
%                 print(f1,'-painter','-depsc', [figuredirectory,fname])
                
                
                %Save everything to data structure
                d.results = results;
                d.fit_plot.x        = x;
                d.fit_plot.pCorrect = fitted_yes;
                d.fit_plot.dprime   = fitted_dprime;
                d.threshold         = threshold; %scaled
                d.slope             = slope; %scaled
                
                output.PYfitdata = d;
                
            end
            
            
            aaa=234;
            
        end
        
        
        
    end
end


%% Plot neurometric data 

for channel = 1:numel(Data.ch)
    for iclu = 1:numel(Data.ch(channel).clu)
        
        if isempty(Data.ch(channel).clu)
            continue
        end
        
        
        
    end
end 



keyboard


% Plot 1) classifier PC fit and d' transformed, 1b) classifier d' formula fit, 2) FR d' formula fit.

% Data.nmdata.classifier
% Data.nmdata.formula

%   ~~ jitters, timebins ~~
% Data.nmdata.classifier.output.PCdata{j}  =  [  stimval  #yesAM  #trials  ]
% Data.nmdata.classifier.output.dprime{j}  =  [ stimval X d' ]
% Data.nmdata.classifier.output.stim{j,1}  =  'jitter'  
% Data.nmdata.classifier.output.stim{j,2}  =  [depths]

% Data.nmdata.formula.FRmat{j}
% Data.nmdata.formula.dpmat{j}






%}







end








