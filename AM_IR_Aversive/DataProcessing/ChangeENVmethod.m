function ChangeENVmethod
% 12/20/19



fn = set_paths_directories();


%%  SUBJECTS

if nargin>0 && exist('select_subject','var')
    if ~iscell(select_subject)
        subjects = {select_subject};
    end
else
    subjects = { 'AAB_265054' 'AAB_265057' 'AAB_265058' 'WWWf_253400' 'WWWlf_253395'}; % nothing good in AAB_265059
end

for subj = 1:numel(subjects)
    
    subject = subjects{subj};
    
    
    %%  SESSIONS
    
    % Get list of sessions to check for sorted data
    
    if nargin>1 && exist('select_session','var')
        if ~iscell(select_session)
            Sessions = {select_session};
        end
    else
        switch subject
            case {'AAB_265054' 'AAB_265058' 'AAB_265057'}
                SpkFns = dir(fullfile(fn.processed,subject,'*AM_TrialData.mat'));
            case {'WWWf_253400' 'WWWlf_253395' }
                SpkFns = dir(fullfile(fn.processed,subject,'*_TrialData.mat'));
        end
        
        Sessions = [];
        for ifn = 1:numel(SpkFns)
            Sessions = [Sessions; extractBetween(SpkFns(ifn).name,'sess-','_TrialData')];
        end
        
    end
    
    
    % Step through each session
    for sess = Sessions'
        
        session = char(sess);
        
        if ~isempty(strfind(session,'-VS'))
            continue
        end
        
        % Load data files
        clear Info TrialData SoundStream SpoutStream Phase0 RateStream
        fprintf('\n***************************************\nLoading %s sess %s...\n',subject,session)
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
        
        
        % Also load epData
        epdatafile = fullfile(fn.raw,subject,Info.epData_fn);
        fprintf(' loading data file %s...',epdatafile)
        clear epData;
        
        epData = load([epdatafile '.mat'],'-mat'); %loads data struct: epData
        
        if ~exist('epData','var')
            keyboard
        elseif isfield(epData,'epData')
            epData = epData.epData;
            disp('(old file) done.')
        end
        disp('done.')
        
        
        % Calculate new SoundStream
        oldSS = SoundStream;
        clear SoundStream SpoutStream Phase0 RateStream
        getPhase0Data;
        
        % Check match between old and new SoundStreams
%         figure; hold on
%         plot(oldSS,'k','LineWidth',2)
%         plot(SoundStream,'LineWidth',2)
        if corr(SoundStream',oldSS')<0.98
            keyboard
        end
        
        
        % Re-save TrialData 
        savename = sprintf('%s_sess-%s_TrialData',SUBJECT,SESS_LABEL);
        save(fullfile(saveDir, savename),'TrialData','SpoutStream','SoundStream','RateStream','Phase0','-v7.3');
        
    end %sess
end %subj


end