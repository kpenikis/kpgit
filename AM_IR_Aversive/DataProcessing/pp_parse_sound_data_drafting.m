function [SoundData,Info] = pp_parse_sound_data_drafting(SoundData,epocs,Info)
%
% called by: pp_processPhys_UMS
% to make SoundData matrix for AM aversive experiment
% 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% IN PROGRESS, AND SWITCHED TO UMS PROGRAM TO EDIT
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Outstanding issues 
%  1) how to easily identify periodic-periodic transitions (iti to rvid
%  when there is a discontinuity/NaN separating them)
    % during raster code. 
    % or possibly: label rate at beginning of each period, 
    % and block type (iti, pdc, ir, unmod) at beginning of each block
    % everything except sound out to timestamps in seconds
%  2) check that silent period is at the beginning
    % during artifact identification?
%  3) what to do if the number of ttls don't line up
%  4) separate this stimulus analysis into a separate function

% find transition from pdc to ir within trial


global fn

% stimfolder = strsplit(Info.stimfiles, '\');
% stimfolder = stimfolder{end};
% stimfolder = fullfile(fn.stim,stimfolder);
% 
% stimfiles = dir(fullfile(stimfolder,'*.mat'));



%Types of sounds
%  silence
%  unmod/Warn
%  periodic iti (4 16 32)
%  periodic half of rV (2 4 8 16 32)
%  irregular half of rV (AC or DB)


% block_id key: 
trial_ID_key = { 'Warn';  '2';   '4';   '8';  '16';  '32';   'AC';  'DB'  };
%   0: silence      1      2      3      4      5      6       7      8

% iti_flag
%  0: trial
%  1: iti

% ir_flag
%  0: periodic
%  1: irregular



% spikes = Spikes.sorted(5);
% spiketimes = round(spikes.spiketimes(spikes.assigns==33) * 1000);  %ms


% Data Table [ n_trials x tags ] 
% each row is a trial, of the types listed above
% columns/tags are:  
   % onset  offset  block_ID  iti_flag  ir_flag  dB  LP  spout_flag  artifact_flag  ... behavior info
TrialData = table;


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

% AT THIS POINT I COULD GET THE TIMESTAMPS ALL IN THE SAME TIMESCALE.
% HOWEVER, I STILL DONT KNOW WHICH TIMESCALE CORRECTLY MATCHES THE
% SPIKETIMES STORED IN THE OTHER CIRCUIT ON THE RZ5.

pdcDur = round(1 * Info.fs_sound); %samples
IRDur  = round(sum(1./[2 4 8 16 32 2 4 8 16 32]) * Info.fs_sound); %samples

nw=0;
nt=0;

for it = 1:numel(trialID)
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % SoundData timescale correction
    % based on ITI onset; assumes epoc timestamps correspond to spiketimes
    
    [~,minidx] = min(abs(1+trialOffset(it) - ITI_up));
    SD_shift = 1+trialOffset(it) - ITI_up(minidx);    % SD idx - SD_shift
    
%     figure;
%     plot(((-1000+trialOnset(it)):(Info.fs_sound+trialOffset(it))),SoundData(2,-SD_shift+((-1000+trialOnset(it)):(Info.fs_sound+trialOffset(it)))),'k')
%     hold on
%     plot(trialOnset(it),0,'c+','MarkerSize',8,'LineWidth',3)
%     plot(trialOffset(it),0,'c+','MarkerSize',8,'LineWidth',3)
%     plot(ITI_up(minidx)+SD_shift,0,'g+','MarkerSize',8,'LineWidth',3)
    
    
    % Which stim make up this trial?
    trialBlocks = strsplit(Info.stimfiles{trialID(it)},'_');
    
    
    %----------------------------------------------------------------------
    % Warn stimulus / unmodulated noise
    if numel(trialBlocks)==2
        
        nw = nw+1;
        
        Warn_Onsets(nw)  = trialOnset(it);
        Warn_Offsets(nw) = trialOffset(it);
        
        Warn_ID(nw)        = 1;
        Warn_ITI_flag(nw) = 0;
        Warn_IR_flag(nw)  = 0;
        
        % Don't have to check stability because params never changed during
        % a Warn trial
        Warn_SPL(nw) = double(SoundData(4,Warn_Onsets(nw)-SD_shift));
        if isfield(epocs,'LPxx')
            Warn_LP(nw) = epocs.LPxx.data(it);
        else
            Warn_LP(nw) = double(SoundData(6,Warn_Onsets(nw)-SD_shift));
        end
        
        % Might want to save time subject left spout
        Warn_Spout_pct(nw) = double(round(100* sum(SoundData(7, (Warn_Onsets(nw):Warn_Offsets(nw))-SD_shift )) / length((Warn_Onsets(nw):Warn_Offsets(nw))))/100);
        
        % Next trial
%         continue
    
    else
    
    
    % Check that duration makes sense
    if abs((trialOnset(it)+pdcDur+IRDur)-trialOffset(it))>1, keyboard, end
    
    nt = nt+1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Periodic block
    
    % onset and offset (samples)
    Pdc_BlockOnsets(nt)  = trialOnset(it);
    Pdc_BlockOffsets(nt) = trialOnset(it) + pdcDur;
    
    % block_ID  iti_flag  ir_flag  dB  LP  spout_flag  artifact_flag  ... behavior info
    Pdc_BlockID(nt) = find(strcmp(trialBlocks{2},trial_ID_key));
    Pdc_ITI_flag(nt) = 0;
    Pdc_IR_flag(nt) = 0;
    
    %%% !! to find corresponding samples for gathering info from SoundData,
    %%%    shifting by the difference between Offset+1 and SD(8,:)/ITI TTL
    %%%    (calculated above)
    
    if numel(unique( SoundData(4, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt))-SD_shift ) ))==1 
        Pdc_SPL(nt) = double(SoundData(4,Pdc_BlockOnsets(nt)-SD_shift));
    else
        Pdc_SPL(nt) = 0;
    end
    
    if isfield(epocs,'LPxx')
        Pdc_LP(nt) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt))-SD_shift ) ))==1 
        Pdc_LP(nt) = double(SoundData(6,Pdc_BlockOnsets(nt)-SD_shift)); 
    else
        Pdc_LP(nt) = 0;
    end
    
    Pdc_Spout_pct(nt) = double(round(100* sum(SoundData(7, (Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt))-SD_shift )) / length((Pdc_BlockOnsets(nt):Pdc_BlockOffsets(nt))))/100);
    
    
    %-  - --  - -  --- - -  -  --    - --  - -- -  - -- --   --- -  -- -- -
    % IR block     
    
    % onset and offset (samples)
    IR_BlockOnsets(nt)  = Pdc_BlockOffsets(nt)+1;
    IR_BlockOffsets(nt) = trialOffset(it);
    
    % block_ID  iti_flag  ir_flag  dB  LP  spout_flag  artifact_flag  ... behavior info
    IR_BlockID(nt) = find(strncmp(trialBlocks{3},trial_ID_key,2));
    IR_ITI_flag(nt) = 0;
    IR_IR_flag(nt) = 1;
    
    %%% !! to find corresponding samples for gathering info from SoundData,
    %%%    shifting by the difference between Offset+1 and SD(8,:)/ITI TTL
    %%%    (calculated above)
    
    if numel(unique( SoundData(4, (IR_BlockOnsets(nt):IR_BlockOffsets(nt))-SD_shift ) ))==1 
        IR_SPL(nt) = double(SoundData(4,IR_BlockOnsets(nt)-SD_shift));
    else
        IR_SPL(nt) = 0;
    end
    
    if isfield(epocs,'LPxx')
        IR_LP(nt) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (IR_BlockOnsets(nt):IR_BlockOffsets(nt))-SD_shift ) ))==1 
        IR_LP(nt) = double(SoundData(6,IR_BlockOnsets(nt)-SD_shift)); 
    else
        IR_LP(nt) = 0;
    end
    
    IR_Spout_pct(nt) = double(round(100* sum(SoundData(7, (IR_BlockOnsets(nt):IR_BlockOffsets(nt))-SD_shift )) / length((IR_BlockOnsets(nt):IR_BlockOffsets(nt))))/100);
    
    
    % To keep in epoc storage timeframe, infer AM phase0s theoretically,
    % from the period durations of the instantaneous rate
    
    
    
    %  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  . 
    % ITI periodic rate  
    
    % Currently not collecting ITIs following Warn trials
    
    % onset and offset (samples)
    ITI_BlockOnsets(nt)  = 1+trialOffset(it);% - ITI_up(minidx)
    if it<length(trialOnset)
        ITI_BlockOffsets(nt) = trialOnset(it+1)-1;
    else
        ITI_BlockOffsets(nt) = size(SoundData,2);
    end
    
    % block_ID  iti_flag  ir_flag  dB  LP  spout_flag  artifact_flag  ... behavior info
    ITI_BlockID(nt)  = find(strncmp( num2str(SoundData(1, ITI_BlockOnsets(nt)+50)), trial_ID_key,2));
    ITI_ITI_flag(nt) = 1;
    ITI_IR_flag(nt)  = 0;
    
    %%% if ITI is longer than the default (ie paused for parameter change,
    %%% detect that and pull clean portion
    %  Leave for now! Only do if there are not enough transition trials
    
    if numel(unique( SoundData(4, (ITI_BlockOnsets(nt):ITI_BlockOffsets(nt))-SD_shift ) ))==1 
        ITI_SPL(nt) = double(SoundData(4,IR_BlockOnsets(nt)-SD_shift));
    else
        ITI_SPL(nt) = 0;
    end
    
    if isfield(epocs,'LPxx')
        ITI_LP(nt) = epocs.LPxx.data(it);
    elseif numel(unique( SoundData(6, (ITI_BlockOnsets(nt):ITI_BlockOffsets(nt))-SD_shift ) ))==1 
        ITI_LP(nt) = double(SoundData(6,ITI_BlockOnsets(nt)-SD_shift)); 
    else
        ITI_LP(nt) = 0;
    end
    
    ITI_Spout_pct(nt) = double(round(100* sum(SoundData(7, (ITI_BlockOnsets(nt):ITI_BlockOffsets(nt))-SD_shift )) / length((ITI_BlockOnsets(nt):ITI_BlockOffsets(nt))))/100);
    
    
    
    end
    
%%
%     %%%%%
%     thisPdcBlock = str2num(trialBlocks{2});
%     thisIRBlock = trialBlocks{3}(1:2);
%     
%     % Get the period offsets corresponding to the periodic block
%     [minval,minidx] = min(abs( Phase0_epoc - trialOnset(it) ));
%     if minval>2, keyboard, end
%     prevPeriodOffset = Phase0_epoc(minidx);
%     prevPeriodRate = epocs.AMrt.data(minidx);
%     
%     % Periodic block onset (samples)
%     PdcBlockOnsets(it) = prevPeriodOffset+1;
%     
%     % Periodic block offset (samples)
%     nPds = pdcDur/1000 * thisPdcBlock;
%     PdcBlockOffsets(it) = Phase0_epoc(minidx + nPds);
%     
%     % Check that duration matches what it should
% %     if abs(ceil(1000*( (PdcBlockOffsets(it) - PdcBlockOnsets(it)) /Info.fs_sound )) - pdcDur) >2
% %         keyboard
% %     end
%     
%     % IR block onset (samples)
%     IRBlockOnsets(it) = PdcBlockOffsets(it)+1;
%     
%     % IR block offset (samples)
%     IRBlockOffsets(it) = Phase0_epoc(minidx + nPds + 10);
%     
%     % Check that duration is reasonable, and ITI onset comes soon after
% %     if abs(ceil(1000*( (IRBlockOffsets(it) - IRBlockOnsets(it)) /Info.fs_sound )) - IRDur) > 5
% %         keyboard
% %     end
%%
    
end

keyboard

% Combine all trial data into Table

% onset  offset  block_ID  iti_flag  ir_flag  dB  LP  spout_flag  artifact_flag  ... behavior info

TrialData.onset   = [Warn_Onsets'; Pdc_BlockOnsets'; IR_BlockOnsets'; ITI_BlockOnsets'];
TrialData.offset  = [Warn_Offsets'; Pdc_BlockOffsets'; IR_BlockOffsets'; ITI_BlockOffsets'];
TrialData.trID    = [Warn_ID'; Pdc_BlockID'; IR_BlockID'; ITI_BlockID'];
TrialData.ITIflag = [Warn_ITI_flag'; Pdc_ITI_flag'; IR_ITI_flag'; ITI_ITI_flag'];
TrialData.IRflag  = [Warn_IR_flag'; Pdc_IR_flag'; IR_IR_flag'; ITI_IR_flag'];
TrialData.SPL     = [Warn_SPL'; Pdc_SPL'; IR_SPL'; ITI_SPL'];
TrialData.LP      = [Warn_LP'; Pdc_LP'; IR_LP'; ITI_LP'];
TrialData.Spout   = [Warn_Spout_pct'; Pdc_Spout_pct'; IR_Spout_pct'; ITI_Spout_pct'];

sortrows(TrialData,'onset')


%%% x why are all the Warn trials at the top (with smallest onsets)?
%%% 1) check collecting params from new saving style session
%%% 2) what things need to be saved as stream:
% phase0 -- could just save a template instead (probably what I was doing
% last time, for Passive experiment)
% spout
% sound waveform?
%%% 3) add behavioral data


return











% Onset of ITI
each_ITI_start = 1+find(diff(SoundData(8,:))==1); %samples


% Full relevant stim data plot
figure;
plot(SoundData(2,:),'k')
hold on
plot((SoundData(1,:)-16)./2./max(SoundData(1,:)),'b')
plot(SoundData(8,:)./2-0.25,'c')
plot(Phase0_epoc,zeros(size(Phase0_epoc)),'r*','MarkerSize',8,'LineWidth',3)
plot((epocs.TTyp.onset)*Info.fs_sound,zeros(size(epocs.TTyp.onset)),'y+','MarkerSize',8,'LineWidth',3)

plot(SoundData(5,:)./2-0.25,'LineWidth',2)
plot(SoundData(6,:)./2-0.25,'LineWidth',2)


[PKS,LOCS,w,p]= findpeaks(double(SoundData(5,:)),'MinPeakProminence',0.9);
plot(LOCS+(w./2),zeros(size(LOCS)),'x','MarkerSize',10,'LineWidth',2)
Phase0_Stream = LOCS+(w./2);


figure;
plot(SoundData(2,:),'k')
hold on
plot(Pdc_BlockOnsets,zeros(size(Pdc_BlockOnsets)),'+','MarkerSize',20,'LineWidth',2)
plot(SoundData(6,:)./2-0.25,'LineWidth',2)

[PKS,LOCS]= findpeaks(double(SoundData(6,:)),'MinPeakProminence',0.9);

LOCS(find(Pdc_BlockOnsets)) - Pdc_BlockOnsets(find(Pdc_BlockOnsets));







% ONSETS OF WARN TRIALS IS A WAY TO ALIGN FOR ALL SESSIONS?
Warn_Onsets = 1+find(diff(SoundData(3,:))==-1);
trialOnset(trialID==12)' - Warn_Onsets(1:end-1);

figure;
plot(SoundData(2,:),'k')
hold on
plot(trialOnset(trialID==12)',-1*ones(size(trialOnset(trialID==12)'))-0.05,'+r','MarkerSize',20,'LineWidth',2)
plot(Warn_Onsets,-1*ones(size(Warn_Onsets))+0.05,'+b','MarkerSize',20,'LineWidth',2)

% plot(spiketimes/1000*Info.fs_sound,-1*ones(size(spiketimes)),'k.')

% linear regression to translate timescales?
p = polyfit(Warn_Onsets(1:length(trialOnset(trialID==12)')),trialOnset(trialID==12)',1)

figure; 
plot(Warn_Onsets(1:length(trialOnset(trialID==12)')),trialOnset(trialID==12)','xk','MarkerSize',10);
hold on
plot(1:1000:length(SoundData),polyval(p,1:1000:length(SoundData)))

% or do same to Phase0 timepoints
p2 = polyfit( [0 Phase0_Stream],[0 Phase0_epoc(1:length(Phase0_Stream))'],1)












% figure;
% plot(SoundData(1,:),'b')
% hold on
% plot(ceil(epocs.AMrt.onset*Info.fs_sound),epocs.AMrt.data,'*r')
% 
% Phase0Changes  = AMpdOffsets(diff(epocs.AMrt.data)~=0);
% SD_rateChanges = find(diff(SoundData(1,:))~=0);
% figure;
% plot(Phase0Changes,'+r')
% hold on
% plot(SD_rateChanges,'+b')
% % There doesn't seem to be a constant shift







% onset of trials
trialOnset  = epocs.rVID.onset'; %MAKE SURE StimTrial GOES HIGH AT PHASE 0 OF FIRST AM PERIOD
% offset of trials
tr_offset_sd = (find(diff(SoundData(8,:))==1)/Info.fs_sound);
% stim ids
trialID = epocs.rVID.data';

% If recording was stopped in the middle of a trial these sizes may not match
if size(trialOnset,2) ~= size(trialOnset,2) 
    keyboard
elseif any((tr_offset_sd-trialOnset)<0) || any((tr_offset_sd-trialOnset)>3)
    keyboard
end


% ITI rates
ITI_rates = unique(SoundData(1,2+find(diff(SoundData(8,:))==1)));
ITI_blockIDs = size(Info.stimfiles,1) + [1 2];
for ir = 1:numel(ITI_rates)
    Info.stimfiles{end+1} = ['ITI_' num2str(ITI_rates(ir))];
end

% SECONDS
figure; 
plot((1:1:size(SoundData,2))/Info.fs_sound,SoundData(8,1:1:end),'r') %iti
hold on
plot((1:1:size(SoundData,2))/Info.fs_sound,SoundData(2,1:1:end)/2+0.5,'Color',[0.6 0.6 0.6]) %sound
plot([trialOnset' tr_offset_sd']',[trialID' trialID']'./max(trialID),'-k') %trial rvids
plot(epocs.AMrt.onset,epocs.AMrt.data./max(epocs.AMrt.data),'g') %AM rate (stored at end of period)
plot((1:1:size(SoundData,2))/Info.fs_sound, SoundData(1,1:1:end)./max(SoundData(1,1:1:end)),'c') %AM rate (sound data)
title('SECONDS')

% SAMPLES
figure; 
plot((1:1:size(SoundData,2)),SoundData(8,1:1:end),'r') %iti
hold on
plot((1:1:size(SoundData,2)),SoundData(2,1:1:end)/2+0.5,'Color',[0.6 0.6 0.6]) %sound
plot([trialOnset' tr_offset_sd']'.*Info.fs_sound, [trialID' trialID']'./max(trialID),'-k') %trial rvids
plot(epocs.AMrt.onset.*Info.fs_sound, epocs.AMrt.data./max(epocs.AMrt.data),'g') %AM rate (stored at end of period)
plot((1:1:size(SoundData,2)), SoundData(1,1:1:end)./max(SoundData(1,1:1:end)),'c') %AM rate (sound data)
title('SAMPLES')



% End goal: add row to SoundData with trial code. rvID plus ITIs, plus 0
% for silence, plus NaN when not an analysis period (e.g. sound cut out, 
% extra period after iti but before next block began)

% Empty data vector
SD_addRow = nan(1,size(SoundData,2));

% Label trial blocks
for ib = 1:numel(trialID)
    ib_samps = ceil([trialOnset(ib) tr_offset_sd(ib)] .*Info.fs_sound);
    SD_addRow(ib_samps(1):ib_samps(end)) = trialID(ib) * ones(1, diff(ib_samps)+1);
end

% Label ITI periods
if ~all(ITI_rates == unique(SoundData(1,2+find(diff(SoundData(8,:))==1))))
    keyboard
end

each_ITI_start = 1+find(diff(SoundData(8,:))==1);
each_ITI_end   = find(diff(SoundData(8,:))==-1);
if numel(each_ITI_start) ~= numel(each_ITI_end), keyboard, end

for ib = 1:numel(each_ITI_start)
    this_ITI = ITI_blockIDs( ITI_rates == double(SoundData(1,1+each_ITI_start(ib))) );
    SD_addRow(each_ITI_start(ib):each_ITI_end(ib)) = this_ITI * ones(1,numel(each_ITI_start(ib):each_ITI_end(ib)));
end
% Samples between end of ITI and start of next block
% round(tr_onset_rv(2:end).*Info.fs_sound - each_ITI_end(1:end-1))

% Label silence (when also drinking)
SD_addRow(1, SoundData(2,:)==0 & SoundData(7,:)==1 ) = 0;

% figure(20);
% plot(SD_addRow,'b')

SoundData = [SoundData; SD_addRow];
Info.sound_rows{end+1} = 'blockID';


% keyboard


% Plot AM rate sources overlaid
figure(21);
plot(SoundData(1,:),'k.')
hold on
plot(epocs.AMrt.onset*Info.fs_sound,epocs.AMrt.data,'b*')

BlockData = nan(size(epocs.AMrt.onset));
for io = 1:numel(epocs.AMrt.onset)
    BlockData(io,1) = SoundData(9,round(epocs.AMrt.onset(io)*Info.fs_sound));
end

plot(epocs.AMrt.onset*Info.fs_sound,BlockData,'r*')


% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


end






