function gatherFullTrialCTTS
% gatherFullTrialCTTS_Speech
%
%   Gathering data for Population analyses...
% 
% KP, 2020-04
%

global fn trN TD Stimuli spkshift 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRIALS
PickTrials   = 'all';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'AC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Unit data files
fn = set_paths_directories;

switch whichStim
    case {'AC' 'DB'}
        TD       = [1500 1000 1000 1000 1000 1000 1937 1937];
        Stimuli  = 1:8;
        trN      = 200;
        UnFN     = 'Units';
    case 'Speech'
        q=load(fullfile(fn.stim,'SpeechStim','TrueSpeechStimDurations'));
        TD       = q.TrueSpeechStimDurations;
        Stimuli  = 1:6;
        trN      = 60;
        UnFN     = 'UnitsVS';
end

q = load(fullfile(fn.processed,UnFN));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

%-------
spkshift = 0; %mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------


%%

Cell_Time_Trial_Stim = []; 
Env_Time_Trial_Stim  = []; 

for iUn = 1:numel(UnitData)
    
    [SpikesTrials,StimTrials] = get_FT_CTTS(UnitData,iUn,PickTrials); 
    
    Cell_Time_Trial_Stim = [Cell_Time_Trial_Stim; SpikesTrials];
    Env_Time_Trial_Stim  = [Env_Time_Trial_Stim; StimTrials];
    
end

% Shorten to max n trials (to save space)
MaxTrs = find(~isnan(mean(mean(mean(Env_Time_Trial_Stim,2,'omitnan'),1,'omitnan'),4,'omitnan')),1,'last');
Cell_Time_Trial_Stim = Cell_Time_Trial_Stim(:,:,1:MaxTrs,:);
Env_Time_Trial_Stim  = Env_Time_Trial_Stim(:,:,1:MaxTrs,:);


% Save data
savedir = fullfile(fn.figs,'PopDynamics');
if ~exist(fullfile(savedir,'Data'),'dir')
    mkdir(fullfile(savedir,'Data'))
end

save(fullfile(savedir,'Data',['CTTS_' whichStim '_' PickTrials]),'Cell_Time_Trial_Stim','Env_Time_Trial_Stim','TD','-v7.3')




%% Modulation power spectrum of speech stimuli

keyboard

fs=1000;

% Get durations of each stimulus
MeanEnv  = permute(mean(mean(Env_Time_Trial_Stim,3,'omitnan'),1,'omitnan'),[4 2 1 3]);
MeanPSTH = permute(mean(mean(Cell_Time_Trial_Stim,3,'omitnan'),1,'omitnan'),[4 2 1 3]);


% Modulation Power Spectrum of all stimuli
CatEnv   = [];
CatPSTH  = [];
for ist =  1:numel(Stimuli) 
    CatEnv  = [CatEnv;  MeanEnv(ist,501:(500+TD(ist)))'];
    CatPSTH = [CatPSTH; MeanPSTH(ist,501:(500+TD(ist)))'];
end

% FFT method #1
[ssf,fws]     = plotMPS(CatEnv,fs);
[ssf_R,fws_R] = plotMPS(CatPSTH,fs);

% FFT method #2
L    = length(CatEnv);
f    = fs*(0:(L/2))'/L;
y_S  = fft(CatEnv);
P2_S = abs(y_S/L);
P1_S = P2_S( 1 : floor(L/2+1) );
P1_S = P1_S.*sqrt(f);
P1_S(2:end-1) = 2*P1_S(2:end-1);

y_R  = fft(CatPSTH);
P2_R = abs(y_R/L);
P1_R = P2_R( 1 : floor(L/2+1) );
P1_R = P1_R.*sqrt(f);
P1_R(2:end-1) = 2*P1_R(2:end-1);


% Plot SIGNAL
figure; 
plot(CatEnv./max(CatEnv),'LineWidth',2)
hold on
plot(CatPSTH./max(CatPSTH),'LineWidth',0.5)
xlim([0 length(CatEnv)])
xlabel('Time (ms)')

print_eps_kp(gcf,fullfile(savedir,'Data',['CatStimPSTH_' whichStim]))


% Plot POWER SPECTRUM
figure; 
xlim([0.25 2^9])
set(gca,'xscale','log')
% plot(ssf,fws./max(fws(ssf>0.25&ssf<64)),'LineWidth',2)
hold on
plot(f,P1_S/max(P1_S),'LineWidth',2)
plot(f,P1_R/max(P1_R),'LineWidth',1)
set(gca,'xtick',[0.5 1 2 4 8 16 32 2^9])
xlabel('Frequency')

print_eps_kp(gcf,fullfile(savedir,'Data',['ModSpectrum_CatStimPSTH_' whichStim]))




end