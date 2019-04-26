function [TrialData,SpoutStream,SoundStream,RateStream,Info] = pp_parse_sound_data(SoundData,epocs,Info)
%
% called by: ProcessPhys_SynKS
% to make TrialData table for AM aversive experiment
% 
% Currently not collecting ITIs following Warn trials
% 
% KP, 2018-03

% Outstanding issues 
%%% add behavioral data
%%% save period dur vectors? or write a simple function to load vector
%%% during analyses


%Types of sounds
%  silence
%  unmod/Warn
%  periodic iti (4 16 32)
%  periodic half of rV (2 4 8 16 32)
%  irregular half of rV (AC or DB)


% stim_id key: 
Info.stim_ID_key = { 'Warn';  '2';   '4';   '8';  '16';  '32';   'AC';  'DB'  };
%   0: silence          1      2      3      4      5      6       7      8



%% Get Stream Data (fs = 1 kHz; 1 ms)

SpoutStream = round( resample(double(SoundData(7,:)),10000,round(Info.fs_sound*10),5) );

% Actually save downsampled RMS 
SoundStream_long = envelope(double(SoundData(2,:)),40,'rms');
SoundStream = resample(SoundStream_long,10000,round(Info.fs_sound*10),5);

RateStream = round(resample(double(SoundData(1,:)),10000,round(Info.fs_sound*10),5));



%%


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% AM period offsets
Phase0_epoc = round(epocs.AMrt.onset * Info.fs_sound); %samples
% [~,LOCS,w,~]= findpeaks(double(SoundData(5,:)),'MinPeakProminence',0.9);
% Phase0_Stream = LOCS+(w./2);
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% Onsets and offsets of trials (StimTrial TTL)
trialOnset   = round(epocs.TTyp.onset * Info.fs_sound); %samples
trialOffset  = round(epocs.TTyp.offset * Info.fs_sound); %samples
trialID      = epocs.rVID.data;
ITI_up       = 1+find(diff(SoundData(8,:))==1);

pdcDur = round(1 * Info.fs_sound); %samples
IRDur  = round(sum(1./[2 4 8 16 32 2 4 8 16 32]) * Info.fs_sound); %samples

nw=0;
nt=0;

for it = 1:numel(trialID)
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % SoundData timescale correction
    % based on ITI onset; assumes epoc timestamps correspond to spiketimes
    
    [~,minidx] = min(abs(1+trialOffset(it) - ITI_up));
%     SD_shift = 1+trialOffset(it) - ITI_up(minidx);    % SD idx - SD_shift
    
%     figure;
%     plot(((-1000+trialOnset(it)):(Info.fs_sound+trialOffset(it))),SoundData(2,-SD_shift+((-1000+trialOnset(it)):(Info.fs_sound+trialOffset(it)))),'k')
%     hold on
%     plot(trialOnset(it),0,'c+','MarkerSize',8,'LineWidth',3)
%     plot(trialOffset(it),0,'c+','MarkerSize',8,'LineWidth',3)
%     plot(ITI_up(minidx)+SD_shift,0,'g+','MarkerSize',8,'LineWidth',3)
    
    
    % Which stim make up this trial?
    blockstring = strsplit(Info.stimfiles{trialID(it)},'_');
    
    
    %----------------------------------------------------------------------
    % Warn stimulus / unmodulated noise
    if numel(blockstring)==2
        
        nw = nw+1;
        
        Warn_Onsets(nw)  = trialOnset(it);
        Warn_Offsets(nw) = trialOffset(it);
        
        Warn_ID(nw)        = 1;
        Warn_ITI_flag(nw) = 0;
        Warn_IR_flag(nw)  = 0;
        
        % Don't have to check stability because params never changed during
        % a Warn trial
        Warn_SPL(nw) = double(SoundData(4,Warn_Onsets(nw)));
        if isfield(epocs,'LPxx')
            Warn_LP(nw) = epocs.LPxx.data(it);
        else
            Warn_LP(nw) = double(SoundData(6,Warn_Onsets(nw)));
        end
        
        % Might want to save time subject left spout
        Warn_Spout_pct(nw) = double(round(100* sum(SoundData(7, (Warn_Onsets(nw):Warn_Offsets(nw)) )) / length((Warn_Onsets(nw):Warn_Offsets(nw))))/100);
        
        
        
        %  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
        % ITI following Warn
                
        % onset and offset (samples)
        if it<length(trialOnset)
            ITI_BlockOffsets(nt+nw) = trialOnset(it+1)-1;
        else
            ITI_BlockOffsets(nt+nw) = size(SoundData,2);
        end
        ITI_BlockOnsets(nt+nw)  = ITI_BlockOffsets(nt+nw) - round(1*Info.fs_sound);
        
%         if ((ITI_BlockOffsets(nt+nw) - ITI_BlockOnsets(nt+nw))/Info.fs_sound) > 1.5
%             keyboard
%         end
        
        ITI_BlockID(nt+nw)  = find(strncmp( num2str(SoundData(1, ITI_BlockOnsets(nt+nw)+50)), Info.stim_ID_key,2));
        ITI_ITI_flag(nt+nw) = 1;
        ITI_IR_flag(nt+nw)  = 0;
        
        if numel(unique( SoundData(4, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) ) ))==1
            ITI_SPL(nt+nw) = double(SoundData(4,ITI_BlockOnsets(nt+nw)+50));
        else
            ITI_SPL(nt+nw) = 0;
        end
        
        if isfield(epocs,'LPxx')
            ITI_LP(nt+nw) = epocs.LPxx.data(it);
        elseif numel(unique( SoundData(6, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) ) ))==1
            ITI_LP(nt+nw) = double(SoundData(6,ITI_BlockOnsets(nt+nw)+50));
        else
            ITI_LP(nt+nw) = 0;
        end
        
        % Check spout only in 1.5 s after onset
        ITI_Spout_pct(nt+nw) = double( round( 100* sum(SoundData(7, (ITI_BlockOnsets(nt+nw) +[0:(1.5*Info.fs_sound)]) )) / (1.5*Info.fs_sound) ) /100);
%         ITI_Spout_pct(nt+nw) = double(round(100* sum(SoundData(7, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) )) / length((ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw))))/100);
        
        
        
        
    else
    
    
    % Check that duration makes sense
    if abs((trialOnset(it)+pdcDur+IRDur)-trialOffset(it))>1, keyboard, end
    
    nt = nt+1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Periodic block
    
    % onset and offset (samples)
    Pdc_BlockOnsets(nt)  = trialOnset(it);
    Pdc_BlockOffsets(nt) = trialOnset(it) + pdcDur;
    
    Pdc_BlockID(nt) = find(strcmp(blockstring{2},Info.stim_ID_key));
    Pdc_ITI_flag(nt) = 0;
    Pdc_IR_flag(nt) = 0;
    
    %%% !! to find corresponding samples for gathering info from SoundData,
    %%%    shifting by the difference between Offset+1 and SD(8,:)/ITI TTL
    %%%    (calculated above)
    
    if numel(unique( SoundData(4, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt)) ) ))==1 
        Pdc_SPL(nt) = double(SoundData(4,Pdc_BlockOnsets(nt)));
    else
        Pdc_SPL(nt) = 0;
    end
    
    if isfield(epocs,'LPxx')
        Pdc_LP(nt) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt)) ) ))==1 
        Pdc_LP(nt) = double(SoundData(6,Pdc_BlockOnsets(nt))); 
    else
        Pdc_LP(nt) = 0;
    end
    
    Pdc_Spout_pct(nt) = double(round(100* sum(SoundData(7, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt)) )) / length((Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt))))/100);
    
    
    %-  - --  - -  --- - -  -  --    - --  - -- -  - -- --   --- -  -- -- -
    % IR block     
    
    % onset and offset (samples)
    IR_BlockOnsets(nt)  = Pdc_BlockOffsets(nt)+1;
    IR_BlockOffsets(nt) = trialOffset(it);
    
    IR_BlockID(nt) = find(strncmp(blockstring{3},Info.stim_ID_key,2));
    IR_ITI_flag(nt) = 0;
    IR_IR_flag(nt) = 1;
    
    if numel(unique( SoundData(4, (IR_BlockOnsets(nt):IR_BlockOffsets(nt)) ) ))==1 
        IR_SPL(nt) = double(SoundData(4,IR_BlockOnsets(nt)));
    else
        IR_SPL(nt) = 0;
    end
    
    if isfield(epocs,'LPxx')
        IR_LP(nt) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (IR_BlockOnsets(nt):IR_BlockOffsets(nt)) ) ))==1 
        IR_LP(nt) = double(SoundData(6,IR_BlockOnsets(nt))); 
    else
        IR_LP(nt) = 0;
    end
    
    IR_Spout_pct(nt) = double(round(100* sum(SoundData(7, (IR_BlockOnsets(nt):IR_BlockOffsets(nt)) )) / length((IR_BlockOnsets(nt):IR_BlockOffsets(nt))))/100);
    
    
    %  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
    % ITI periodic rate  
        
    % onset and offset (samples)
    ITI_BlockOnsets(nt+nw)  = 1+trialOffset(it);% - ITI_up(minidx)
    if it<length(trialOnset)
        ITI_BlockOffsets(nt+nw) = trialOnset(it+1)-1;
    else  %if it's the last stim in the recording session
        ITI_BlockOffsets(nt+nw) = size(SoundData,2);
    end
    
    %%% If ITI is longer than the default (ie paused for parameter change,
    %%% detect that and pull clean portion
    %%% Leave as is for now! Only do if there are not enough transition trials
%     if ((ITI_BlockOffsets(nt+nw) - ITI_BlockOnsets(nt+nw))/Info.fs_sound) > 1.5
%         keyboard
%     end
    
    ITI_BlockID(nt+nw)  = find(strncmp( num2str(SoundData(1, ITI_BlockOnsets(nt+nw)+50)), Info.stim_ID_key,2));
    ITI_ITI_flag(nt+nw) = 1;
    ITI_IR_flag(nt+nw)  = 0;
    
    if numel(unique( SoundData(4, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) ) ))==1 
        ITI_SPL(nt+nw) = double(SoundData(4,ITI_BlockOnsets(nt+nw)+50));
    else
        ITI_SPL(nt+nw) = 0;
    end
    
    if isfield(epocs,'LPxx')
        ITI_LP(nt+nw) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) ) ))==1 
        ITI_LP(nt+nw) = double(SoundData(6,ITI_BlockOnsets(nt+nw)+50)); 
    else
        ITI_LP(nt+nw) = 0;
    end
    
    % Check spout only in 1.5 s after onset
    ITI_Spout_pct(nt+nw) = double( round( 100* sum(SoundData(7, (ITI_BlockOnsets(nt+nw) +[0:(1.5*Info.fs_sound)]) )) / (1.5*Info.fs_sound) ) /100);
%     ITI_Spout_pct(nt+nw) = double( round( 100* sum(SoundData(7, (ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw)) )) / length((ITI_BlockOnsets(nt+nw):ITI_BlockOffsets(nt+nw))))/100);
    
    end %if numel(blockstring)==2  ->  Warn stim or other 
    
end %it



% Get timestamps of Silent period at beginning of session (for baseline)

win = round(1*Info.fs_sound);
coeffFilt   = ones(1, win)/win;
fDelay    = (length(coeffFilt)-1)/2;

Silence_Onset  = round(Info.fs_sound) + find(abs(filter(coeffFilt, 1, double(SoundData(7,:)))-1)<1e-3,1,'first')-fDelay; %when spout goes high: SoundData(7,:) + 1 second
Silence_Offset = find(SoundData(2,:),1,'first'); %when sound begins - SoundData(2,:)
Silence_Spout_pct = double(round(100* sum(SoundData(7, (Silence_Onset:Silence_Offset) )) / length((Silence_Onset:Silence_Offset)))/100);



%% Make TrialData
% Combine all trial data into Table (and convert samples to ms)

% Data Table [ n_trials x tags ] 
% each row is a trial, of the types listed above
% columns/tags are:  
% onset(ms)  offset(ms)  trial_ID  iti_flag  ir_flag  dB  LP  spout_flag ... behavior info ... artifact_flag

TrialData = table;
TrialData.onset   = round([Silence_Onset;  Warn_Onsets';  Pdc_BlockOnsets';  IR_BlockOnsets';  ITI_BlockOnsets'] ./Info.fs_sound.*1000); %ms
TrialData.offset  = round([Silence_Offset; Warn_Offsets'; Pdc_BlockOffsets'; IR_BlockOffsets'; ITI_BlockOffsets']./Info.fs_sound.*1000); %ms
TrialData.trID    = [0; Warn_ID'; Pdc_BlockID'; IR_BlockID'; ITI_BlockID'];
TrialData.ITIflag = [0; Warn_ITI_flag'; Pdc_ITI_flag'; IR_ITI_flag'; ITI_ITI_flag'];
TrialData.IRflag  = [0; Warn_IR_flag'; Pdc_IR_flag'; IR_IR_flag'; ITI_IR_flag'];
TrialData.SPL     = [0; Warn_SPL'; Pdc_SPL'; IR_SPL'; ITI_SPL'];
TrialData.LP      = [0; Warn_LP'; Pdc_LP'; IR_LP'; ITI_LP'];
TrialData.Spout   = [Silence_Spout_pct; Warn_Spout_pct'; Pdc_Spout_pct'; IR_Spout_pct'; ITI_Spout_pct'];

TrialData = sortrows(TrialData,'onset');


%% Update stimID key for labeling in analyses

Info.stim_ID_key = { 'Warn';  '2Hz';   '4Hz';   '8Hz';  '16Hz';  '32Hz';   'AC';  'DB'  };
%   0: silence          1       2        3        4       5        6        7       8



end






