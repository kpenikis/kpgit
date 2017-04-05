function Data = getALL_nmData(subj,sess)


% Set some basic variables 
clear  fn subject session channel clu paramdir
global fn subject session channel clu paramdir

subject = subj; 
session = sess;
fn = set_paths_directories(subject,session);

set(0,'DefaultTextInterpreter','none')

% Add psignifit path and turn off warnings
addpath('/Users/kpenikis/Documents/MATLAB/psignifit_git')
warning('off','stats:nlinfit:ModelConstantWRTParam')
warning('off','stats:nlinfit:IllConditionedJacobian')
warning('off','stats:nlinfit:Overparameterized')
warning('off','MATLAB:rankDeficientMatrix')
warning('off','psignifit:ThresholdPCchanged');

% Load Data file for session
dataname = sprintf( '%s_sess-%s_Data',subject,session);
if nargin<3
    load(fullfile(fn.processed,subject,dataname))
end

%%  Calculate neurometric data

% Go through each cluster and call neurometric programs
allclusters = fieldnames(Data);

for unit = allclusters'
    
    clear cludata; clu_nmdata = struct;
    channel = Data.(unit{:}).labels(1,3);
    clu     = Data.(unit{:}).labels(1,1);
    Stimuli = Data.(unit{:}).stimdata;
    
    close all
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~   ch %i, clu %i   ~~~~~~~~~~~~~~~~~~~~~~~~\n\n',channel,clu)
        
    % For each stimulus set presented in this session, and for each set of
    % acoustic parameters, get neurometric data.
    for is = 1:numel(Stimuli)
        
        for ip = 1:numel(Stimuli(is).pars)
            
            fprintf(' ---      pars (HP LP dB AMrate):   %s      --- \n\n',num2str(Stimuli(is).pars(ip).pars))
            
            % Set folder name for saving figures
            paramdir  = sprintf('AM%iHz_%idB_%i-%i_blk%i',...
                Stimuli(is).pars(ip).pars(1,4),...
                Stimuli(is).pars(ip).pars(1,3),...
                Stimuli(is).pars(ip).pars(1,1),...
                Stimuli(is).pars(ip).pars(1,2),...
                Stimuli(is).block);
            
            
            %% Run discrimination analyses 
            %  depending on what feature can be discriminated in this stim set
            
            %%%%%%%%%%%%%%
            %%% JITTER discrimination only
            if isfield(Stimuli(is).pars(ip),'nm_jitter') && ~isfield(Stimuli(is).pars(ip),'nm_depth')
                                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,clu)
                nm_jitter = collect_nmData(Stimuli(is).pars(ip).nm_jitter,'jitter');
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to this cluster's nm data struct
                clu_nmdata(is).pars(ip).nm_jitter = nm_jitter;
                
                
            
            %%%%%%%%%%%%%%
            %%% DEPTH discrimination only
            elseif isfield(Stimuli(is).pars(ip),'nm_depth') && ~isfield(Stimuli(is).pars(ip),'nm_jitter')
                
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,clu)
                nm_depth = collect_nmData(Stimuli(is).pars(ip).nm_depth,'depth');
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to this cluster's nm data struct
                clu_nmdata(is).pars(ip).nm_depth = nm_depth;
                
                                
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% JITTER and DEPTH discrimination
            elseif isfield(Stimuli(is).pars(ip),'nm_jitter') && isfield(Stimuli(is).pars(ip),'nm_depth')
                
                keyboard
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,clu)
                nm_depth = collect_nmData(Stimuli(is).pars(ip).nm_depth,'depth');
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to this cluster's nm data struct
                clu_nmdata(is).pars(ip).nm_depth = nm_depth;
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,clu)
                nm_jitter = collect_nmData(Stimuli(is).pars(ip).nm_jitter,'jitter');
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to this cluster's nm data struct
                clu_nmdata(is).pars(ip).nm_jitter = nm_jitter;
                
                
                
                
            end
        end
    end
    
    % * * * * * * * * * * * * * * * * *
    % Save cluster structure to file
    % * * * * * * * * * * * * * * * * *
    disp('## Saving file for this cluster.')
    eval(sprintf('%s = clu_nmdata;',unit{:}));
    cluefilename = sprintf( '%s_sess-%s_%s',subject,session,unit{:});
    save(fullfile(fn.sess_data,cluefilename),unit{:},'-v7.3');
    
    % Clear from workspace
    eval(sprintf('clear(''%s'');',unit{:}));
    clear clu_nmdata
    
end


%%


end












