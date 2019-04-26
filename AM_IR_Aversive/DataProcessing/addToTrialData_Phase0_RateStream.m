function addToTrialData_Phase0_RateStream
%
%  save_Phase0_RateStream
%
%  KP, 2019-04
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

for ii = 1:size(SESSIONS,1)
    
    subject   = SESSIONS{ii,1}{:};
    session   = SESSIONS{ii,2}{:};
    
    % Load Info and TrialData files
    clear Info TrialData SoundStream SpoutStream RateStream Phase0 
    filename = sprintf( '%s_sess-%s_Info'      ,subject,session); load(fullfile(fn.processed,subject,filename));
    filename = sprintf( '%s_sess-%s_TrialData' ,subject,session); load(fullfile(fn.processed,subject,filename));
%     clear RateStream Phase0
    if exist('RateStream','var') && exist('Phase0','var')
        continue
    end
    
    % Load epData file
    fprintf('Loading epData for %s %s... \n',subject,session)
    try
        load(fullfile(fn.raw,subject,Info.epData_fn));
    catch
        try
            load(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '.mat']));
        catch 
            keyboard
        end
    end
    
    
    % Get RateStream (fs = 1 kHz; 1 ms)
    idxms = round( 1 : Info.fs_sound/1000 : length(epData.streams.rVrt.data(1,:)) );
    RateStream = double(epData.streams.rVrt.data(1,idxms));
    
    
    % Get Phase0 timestamps (fs = 1 kHz; 1 ms precision)
    [~,AMMIN] = findpeaks(-SoundStream./max(SoundStream),'MinPeakProminence',0.5); %ms of RMS minima
    Phase0 = [AMMIN; RateStream(AMMIN+10)];
    
    % Remove outlier periods
    pctDiff = abs( (RateStream(AMMIN(1:end-1)+10) - round(1000./diff(AMMIN)) ) ./ RateStream(AMMIN(1:end-1)+10) );
    
    for irate=[2 4 8 16 32]
        idx = find(Phase0(2,:)==irate);
        fprintf('   %i of %i %i Hz pds removed\n', sum(pctDiff(idx(1:end-1))>0.1),numel(idx),irate)
    end
    
    Phase0(2,find(pctDiff>0.1)) = nan;
    
%     for irate=[2 4 8 16 32]
%         hf(irate) = figure; hold on
%         idx = find(Phase0(2,:)==irate);
%         for imin = idx
%             if imin==length(Phase0), continue, end
%             plot(SoundStream(Phase0(1,imin):Phase0(1,imin+1)),'k','LineWidth',2)
%             xlim([1 10+1000/irate])
%             ylim([0 1.5])
%         end
%     end
%     for irate=[2 4 8 16 32]
%         figure(hf(irate)); hold on
%         idx = find(Phase0(2,:)==irate);
%         for imin = idx
%             if imin==length(Phase0), continue, end
%             
%             if pctDiff(imin)>0.1
%                 plot(SoundStream(Phase0(1,imin):Phase0(1,imin+1)),'r','LineWidth',2)
%             end
%         end
%     end
    
    
    try 
        
    % Re-save TrialData
    save(fullfile(fn.processed,subject,filename),'TrialData','SpoutStream','SoundStream','RateStream','Phase0','-v7.3');
    
    % Re-save epData with struct format
    fprintf('                            saving... ')
    try
        save(fullfile(fn.raw,subject,Info.epData_fn),'-struct','epData','-v7.3')
    catch
        save(fullfile(fn.raw,subject,['Block-' num2str(Info.block) '.mat']),'-struct','epData','-v7.3')
    end
    fprintf('done.\n')
    
    catch
        keyboard
    end
end





end




