function Data = ap_neurometricALL(subj,sess,Data)


global subject session channel 
subject = subj; session = sess;

%~~~~~~~~~~~~~~~
SKIPCLASS = 1;
%~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~
SKIPFIT   = 1;
%~~~~~~~~~~~~~~~

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

addpath(genpath('analysis_modules'),'helpers')
addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')

warning('off','stats:nlinfit:ModelConstantWRTParam')
warning('off','stats:nlinfit:IllConditionedJacobian')
warning('off','stats:nlinfit:Overparameterized')
warning('off','MATLAB:rankDeficientMatrix')
warning('off','psignifit:ThresholdPCchanged');

% Load Data file for session
datadir  = '/Users/kpenikis/Documents/SanesLab/Data/AMJitter/ProcessedData';
filename = sprintf( '%s_sess-%s_Data',subject,session); 
if nargin<3
load(fullfile(datadir,subject,filename))
end

%%  Calculate neurometric data 


% Step through each clu of each channel
for channel = 11%:numel(Data.ch)
    
    % Remove the field to start over clean
    if isfield(Data.ch(channel).clu,'nmdata') 
        Data.ch(channel).clu = rmfield([Data.ch(channel).clu],'nmdata');
    end
    
    for iclu = 1:numel(Data.ch(channel).clu)
        
        if isempty(Data.ch(channel).clu)
            continue
        end
        
        cludata = Data.ch(channel).clu(iclu);
        
        close all
        fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~   ch %i, clu %i   ~~~~~~~~~~~~~~~~~~~~~~~~\n\n',channel,Data.ch(channel).clu(iclu).label(1))
        
        %  Calculate neurometric data if not yet run, or if rerun specified.
            
            % For each block of experiments, determine what parameter
            % was varied and can be discriminated.
            for is = 1:numel(cludata.analyses.stim)
                for ip = 1:numel(cludata.analyses.stim(is).pars)
                    
                    fprintf(' ---      pars (HP LP dB AMrate):   %s      --- \n\n',num2str(cludata.analyses.stim(is).pars(ip).pars))
                    
                    
                    jitters  = cellfun(@length,{cludata.analyses.stim(is).pars(ip).depths.jitter},'UniformOutput',false);
                    njitters = vertcat(jitters{:});
                    
                    %%% JITTER discrimination only
                    if( length(njitters)<3 && any(njitters>=3) )
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,Data.ch(channel).clu(iclu).label(1))
                        [nm_jitter,cludata] = collect_nmData(cludata,is,ip,'jitter');
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        % Add to Data struct and save file
                        cludata.analyses.stim(is).pars(ip).nm_jitter = nm_jitter;
                        Data.ch(channel).clu(iclu) = cludata;
                        
                        
                    %%% DEPTH discrimination only
                    elseif( length(njitters)>=3 && any(njitters<3) )
                        
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,Data.ch(channel).clu(iclu).label(1))
                        [nm_depth,cludata] = collect_nmData(cludata,is,ip,'depth');
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        % Add to Data struct and save file
                        cludata.analyses.stim(is).pars(ip).nm_depth = nm_depth;
                        Data.ch(channel).clu(iclu) = cludata;
                        

                    %%% JITTER and DEPTH discrimination    
                    elseif( length(njitters)>=3 && any(njitters>=3) )
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,Data.ch(channel).clu(iclu).label(1))
                        [nm_depth,cludata] = collect_nmData(cludata,is,ip,'depth');
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        % Add to Data struct and save file
                        cludata.analyses.stim(is).pars(ip).nm_depth = nm_depth;
                        Data.ch(channel).clu(iclu) = cludata;                        
                        
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,Data.ch(channel).clu(iclu).label(1))
                        [nm_jitter,cludata] = collect_nmData(cludata,is,ip,'jitter');
                        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        
                        % Add to Data struct and save file
                        cludata.analyses.stim(is).pars(ip).nm_jitter = nm_jitter;
                        Data.ch(channel).clu(iclu) = cludata;
                        
                        
                    end
                end
            end
            
            % Save Data structure
            save(fullfile(datadir,subject,filename),'Data','-v7.3');
            
    end
end

end











