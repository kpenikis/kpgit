function addToTrialData_Phase0_RateStream
%
%  save_Phase0_RateStream
%
%  KP, 2019-04, updated 2019-06 
%

%%
        %     fs = epData.streams.rVrt.fs (1/5 of phys data)
        % epData.streams.rVrt.data
        %  (1,:) = Instantaneous AM rate  (of the period that is just ending)
        %  (2,:) = Sound output          
        %  (3,:) = AM depth
        %  (4,:) = dB SPL
        %  (5,:) = HP
        %  (6,:) = LP
        %  (7,:) = Spout TTL
        %  (8,:) = ITI TTL
        
        % falling edge of InTrial:
        %   epData.epocs.RCod  =  response code
        % StimTrial TTL:
        %   epData.epocs.TTyp  =  TrialType (0 safe, 1 warn)
        %   epData.epocs.Opto  =  optostim
        %   epData.epocs.rVID  =  rateVec_ID
        % rising edge of PHASE0:
        %   epData.epocs.AMrt  =  AMRATE
           
%%

global fn
close all

fn = set_paths_directories('','',1);

% Load Unit data files
q = load(fullfile(fn.processed,'Units'));
UnitInfo = q.UnitInfo;
clear q

SESSIONS = unique(UnitInfo(:,1:2));
% ISSUE WITH Jan31-AM  -- CORRUPTED

for ii = 1:size(SESSIONS,1)
    
    subject   = SESSIONS{ii,1}{:};
    session   = SESSIONS{ii,2}{:};
    
    % Load Info and TrialData files
    clear Info TrialData
    filename = sprintf( '%s_sess-%s_Info'      ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); load(fullfile(fn.processed,subject,filename));
    
    % Change this depending on how you want to update the stim data
    clear SoundData SoundStream SpoutStream RateStream Phase0 
    
    % Load epData file
    clear epData
    fprintf('Loading epData for %s %s... \n',subject,session)
    try
        epData = load([fullfile(fn.raw,subject,Info.epData_fn) '.mat'],'-mat');
    catch
        try
            epData = load(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '.mat']),'-mat');
        catch 
            keyboard
        end
    end
    if ~exist('epData','var')
        fprintf('COULDN''T LOAD epData..skipping..\n')
        continue
    elseif isfield(epData,'epData')
        epData = epData.epData;
        disp('(old file) done.')
    end
    disp('done.')
    
    
    %%%%%%%%%%%%%%
    getPhase0Data;
    %%%%%%%%%%%%%%
    
    try 
        
    % Re-save TrialData
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); 
    save(fullfile(fn.processed,subject,filename),'TrialData','SpoutStream','SoundStream','RateStream','Phase0','-v7.3');
    % Re-save Info  as well
    filename = sprintf( '%s_sess-%s_Info',subject,session);
    save(fullfile(fn.processed,subject,filename),'Info','-v7.3');
    
    % Re-save epData with struct format
%     fprintf('                            saving... ')
%     try
%         save(fullfile(fn.raw,subject,Info.epData_fn),'-struct','epData','-v7.3')
%     catch
%         save(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '.mat']),'-struct','epData','-v7.3')
%     end
%     fprintf('done.\n')
%     
    catch
        keyboard
    end
end





end




