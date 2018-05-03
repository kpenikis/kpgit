function coldcall_identify_artifact(subject,session,Info,Phys,TrialData)
% coldcall_identify_artifact(subject,session,Info,Phys,TrialData)
%  Calls identify_artifact_trials, for data that has been run through
%  pp_processPhys. 
%  Required inputs are subject and session. If Info, Phys, TrialData are
%  already in workspace, can include them as input vars to save time
%  loading. 
%  KP, 2018-02
%

fn = set_paths_directories(subject,session,1);
if nargin<5
    fprintf('loading data...\n')
    filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_Phys'     ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
end

Info = identify_artifact_trials(Info,Phys,TrialData);


%% Save data

fprintf('\nsaving data...\n')

try
    savedir = fullfile(fn.processed,subject);
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    
    % Save Info structure
    savename = sprintf('%s_sess-%s_Info',subject,session);
    save(fullfile( savedir, savename),'Info','-v7.3');
    
catch
    keyboard
end


end