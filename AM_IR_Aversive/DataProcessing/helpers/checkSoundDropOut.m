

fprintf(' checking for sound drop out...\n')

data = double(SoundData(2,:)); 
 
data_rms = envelope(data,round(Info.fs_sound*2),'rms');
 
samp_beg = Info.fs_sound*TrialData.onset(2)/1000;
samp_end = Info.fs_sound*TrialData.onset(end)/1000;
 
thresh = std(data_rms)*4 + mean(data_rms);
 
flagged = find(data_rms>thresh);
if ~isempty(flagged)
    buffer = round(Info.fs_sound*0.5);
    flagged = (flagged(1)-buffer):(flagged(end)+buffer);
end

% figure; plot(envelope(data,round(Info.fs_sound*2),'rms'))
% hold on
% plot(flagged,thresh*ones(size(flagged)),'r')

data(flagged)    = 0; % can't be nans for envelope function

% Remake sound stream
SoundStream_long = envelope(data,40,'rms');
SoundStream      = resample(SoundStream_long,10000,round(Info.fs_sound*10),5);

% Also resample sound flag to mark artifact trials later
flag_01          = zeros(size(data));
flag_01(flagged) = 1;
SoundFlag        = resample(flag_01,10000,round(Info.fs_sound*10),5);


%% Flag the trials that have no sound

ArtifactFlags = zeros(size(TrialData,1),1);

for it = 1:size(TrialData,1)
    if sum( ismember(TrialData.onset(it):TrialData.offset(it), find(SoundFlag)) )>0 && it>1
        ArtifactFlags(it,1) = 1;
    end
end

%%Save the artifact flags into the Info struct
Info.soundflag.trials = find(ArtifactFlags(:,1));








 