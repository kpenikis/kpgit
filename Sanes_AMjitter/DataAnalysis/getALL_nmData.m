function Data = getALL_nmData(subj,sess)


clear  fn subject session channel paramdir
global fn subject session channel paramdir


subject = subj; 
session = sess;
fn = set_paths_directories(subject,session);

%~~~~~~~~~~~~~~~
SKIPCLASS = 0;
%~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~
RERUNALL  = 1;
%~~~~~~~~~~~~~~~

set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextInterpreter','none')

% Add psignifit path and turn off some warnings
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
    
    if ~strncmpi(unit{:},'ch',2), continue, end
    
    addJitter = [];  addDepth  = [];
    
    % Load cluster data file
    clear clustruct Data_clu
    cluefilename = sprintf( '%s_sess-%s_%s',subject,session,unit{:});
    clustruct = load(fullfile(fn.sess_data,cluefilename));
    clustruct = clustruct.(unit{:});
    channel = clustruct.labels(1,3);
    Data_clu = Data.(unit{:});
    
    close all
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~   ch %i, clu %i   ~~~~~~~~~~~~~~~~~~~~~~~~\n\n',channel,clustruct.labels(1,1))
        
    % For each stimulus set presented in this session, and for each set of
    % acoustic parameters, get neurometric data.
    for is = 1:numel(clustruct.block)
        
        for ip = 1:numel(clustruct.block(is).pars)
            
            % Remove the relevant fields to start over clean
            if RERUNALL
                if      isfield(clustruct.block(is).pars(ip).stimvals,'raster_idxs')
                        clustruct.block(is).pars(ip).stimvals = rmfield(clustruct.block(is).pars(ip).stimvals,'raster_idxs');
                end
%                 if      isfield(clustruct.block(is).pars,'nm_jitter')
%                         clustruct.block(is).pars = rmfield(clustruct.block(is).pars,'nm_jitter');
%                         
%                 elseif  isfield(clustruct.block(is).pars,'nm_depth')
%                         clustruct.block(is).pars = rmfield(clustruct.block(is).pars,'nm_depth');
%                 end
            end
            
            fprintf(' ---      pars (HP LP dB AMrate):   %s      --- \n\n',num2str(clustruct.block(is).pars(ip).pars))
            
            % Set folder name for saving figures
            bk_str = num2str(unique([clustruct.block(is).block]));
            bk_str = bk_str(~isspace(bk_str));
            paramdir  = sprintf('AM%iHz_%idB_%i-%i_blk%s',...
                clustruct.block(is).pars(ip).pars(1,4),...
                clustruct.block(is).pars(ip).pars(1,3),...
                clustruct.block(is).pars(ip).pars(1,1),...
                clustruct.block(is).pars(ip).pars(1,2), bk_str);
            
            
            %% Get the number of jitters and depths in this stimulus set
            njitters  = cell2mat(cellfun(@length,{clustruct.block(is).pars(ip).stimvals.jitter},'UniformOutput',false));
            
            %%%%%%%%%%%%%%
            %%% JITTER discrimination only
            if( length(njitters)<3 && any(njitters>=3) )
                                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,clustruct.labels(1,1))
                [nm_jitter,clustruct,Data_clu] = collect_nmData(clustruct,is,ip,'jitter',Data_clu);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to Data struct and save file
                clustruct.block(is).pars(ip).nm_jitter = nm_jitter;
                
                addJitter = [addJitter [is;ip]];
                
            
            %%%%%%%%%%%%%%
            %%% DEPTH discrimination only
            elseif( length(njitters)>=3 && any(njitters<3) )
                
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,clustruct.labels(1,1))
                [nm_depth,clustruct,Data_clu] = collect_nmData(clustruct,is,ip,'depth',Data_clu);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to Data struct and save file
                clustruct.block(is).pars(ip).nm_depth = nm_depth;
                
                addDepth  = [addDepth [is;ip]];
                
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% JITTER and DEPTH discrimination
            elseif( length(njitters)>=3 && any(njitters>=3) )
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric DEPTH data for channel %i clu %i\n',channel,clustruct.labels(1,1))
                [nm_depth,clustruct,Data_clu] = collect_nmData(clustruct,is,ip,'depth',Data_clu);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to Data struct and save file
                clustruct.block(is).pars(ip).nm_depth = nm_depth;
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                fprintf('Getting neurometric JITTER data for channel %i clu %i\n',channel,clustruct.labels(1,1))
                [nm_jitter,clustruct,Data_clu] = collect_nmData(clustruct,is,ip,'jitter',Data_clu);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                % Add to Data struct and save file
                clustruct.block(is).pars(ip).nm_jitter = nm_jitter;
                
                addJitter = [addJitter [is;ip]];
                addDepth  = [addDepth  [is;ip]];
                
                
            end
        end
    end
    

    % Save cluster structure to file
    disp('## Saving file for this cluster.')
    eval(sprintf('%s = clustruct;',unit{:}));
    save(fullfile(fn.sess_data,cluefilename),unit{:},'-v7.3');
    
    % Check out Data clu also -- particularly if raster was added
    Data.(unit{:}) = Data_clu;
    
    eval(sprintf('clear(''%s'');',unit{:}));
    
end


%%

% Add marker for stim type
if ~isempty(addJitter)
    Data.jitter = addJitter;
end
if ~isempty(addDepth)
    Data.depth  = addDepth;
end

% Save Data file
save(fullfile(fn.processed,subject,dataname),'Data','-v7.3');



end












