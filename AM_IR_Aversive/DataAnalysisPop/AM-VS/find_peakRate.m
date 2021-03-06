function [allTS, env, peakRate, peakEnv, minRatio] = find_peakRate(sound, soundfs, onsOff, envtype)
% function [env, peakRate, peakEnv] = find_peakRate(sound, soundfs, onsOff, envtype)
% Inputs: 
%   sound - time x 1, sound waveform
%   soundfs - 1x1, sampling frequency of sound
%   onsOff  - 1 x 2 times of stimulus onset and offset in sound (in seconds)
%   envtype: 'loudness' (default), or 'broadband': specific loudness envelope or broadband envelope
% Output: 
%   env - amplitude envelope of input
%   peakRate - discrete time series of peakRate events in envelope
%   peakEnv - discrete time series of peakEnv events

% this function creates the speech envelope and a timeseries of discrete
% peakRate events, as used in Oganian & Chang, 2019, A speech envelope
% landmark for syllable encoding in human superior temporal gyrus, Science
% Advances, xxxx

% (c) Yulia Oganian, Oct 2019
% yulia.oganian@ucsf.edu

%% initialize
% need to know when signal contains speech to remode landmark events
% outside speech. this is t
if nargin<3
    onsOff= [1/soundfs length(sound)/soundfs];
end

if nargin<4, envtype = 'rms'; end
% if nargin<4, envtype = 'broadband'; end

% if set to 1, landmark series will be cleaned up to contain only a single
% peakRate event in each envelope cycle, defined as envelope
% trough-to-trough
cleanup_flag = 0;

minRatio     = 1;

%% get envelope

envfs    = 1000;
filtprop = envfs/soundfs/20;
mpd      = 1 / soundfs * envfs; %minpeakdistance: convert to env fs
% mpd      = 6; %minpeakprominence !

switch envtype
    case 'loudness' %% specific loudness
        env = Syl_1984_Schotola(sound, soundfs);        
    case 'broadband' %% broadband envelope        
        rectsound  = abs(sound);
        [b,a] = butter(4, filtprop); %10/(soundfs/2)
        cenv = filtfilt(b, a, rectsound);
        downsenv = resample(cenv, (1:length(cenv))/soundfs, envfs);
        downsenv(downsenv <0) =0;
        env = downsenv;
    case 'rms'
        env = sound;
        envfs = soundfs;
end

%% get landmarks in envelope
%
% vectors marking speech onset and offset times.
onsOffenv = zeros(2,length(env));
onsOffenv(1,ceil(onsOff(1)*envfs))  = 1;
onsOffenv(2,round(onsOff(2)*envfs)) = 1;


allTS = find_landmarks(env, onsOffenv, cleanup_flag, mpd, minRatio); 
peakEnv  = allTS(4,:);
peakRate = allTS(6,:);

% figure; hold on
% plot(env)
% plot(peakRate)
% plot(diff(env))

% figure; 
% histogram(peakRate,logspace(-3,-1,100))
% ylim([0 15])
% xlim([1e-3 1e-1])
% set(gca,'xscale','log')

end

%% specific loudness function for envelope extraction
function Nm = Syl_1984_Schotola(p,fs)
%
%This identifies vowell nuclei according to the method of Schotola (1984).
%The calculation of modifed total loudness is identical to that of
%Syl_1979_Zwicker_etal_v1.m (e.g., Zwicker et al. 1979), but some
%additional rules/improvements are made here. This is continuation of work
%by Ruske and Schotola (1978, 1981), but this is the best-specified.
%
%INPUTS:
%p [N x 1]: speech waveform, p(t)
%fs [num]: sample rate (Hz)
%
%OUTPUTS:
%Nm [N x 1]: modified total loudness, Nm(t), used to identify VN and B


% written by Eric Edwards 
% adopted by Yulia Oganian, yulia.oganian@ucsf.edu

p = p(:);
N = length(p);
tN = (0:N-1)'./fs;
T = 1/fs;

%Loudness functions will be sampled at 100 Hz
sr = 100;
N1 = fix(N*sr/fs);
t1 = (0:N1-1)'./sr;

%Critical-band filters are applied in freq-domain
p = fft(p,N);
frqs = (0:N/2)*(fs/N); %FFT positive freqs
nfrqs = length(frqs);

%Set up for cricial-band filter bank (in freq-domain)
z = 13*atan(.76*frqs/1000) + 3.5*atan(frqs/7500).^2; %Bark (Zwicker & Terhardt 1980)
%CB=25+75*(1+1.4*(frqs/1000).^2).^.69; %Critical bandwidth (Zwicker & Terhardt 1980)
z = z(2:nfrqs-1); %to set DC and Nyquist to 0

%Set up for simple (RC) smoothing with 1.3-ms time-constant
tau = 0.0013; r = exp(-T/tau);
b1 = 1-r; a1 = poly(r);

F = zeros([N 1],'double'); %will hold critical-band filter shape
czs = 1:22; nczs = length(czs);
Nv = zeros ([N1 nczs],'double'); %will hold the specific loudnesses
for ncz=1:nczs, cz = czs(ncz);
    F(2:nfrqs-1) = 10.^(.7-.75*((z-cz)-.215)-1.75*(0.196+((z-cz)-.215).^2));
    %F = F./sum(F);
    Lev = real(ifft(p.*F,N));
    Lev = Lev.^2; %square-law rectification (Vogel 1975)
    Lev = filter(b1,a1,Lev); %smoothing (1.3-ms time-constant)
    Lev = flipud(filter(b1,a1,flipud(Lev)));
    %the last line makes zero-phase smoothing, comment out to leave causal
    
    Lev = log(Lev); %logarithmic nonlinearity; this is now "excitation level"
    %The "specific loudness", Nv(t), is made by "delog" of .5*Lev(t)
    Nv(:,ncz) = interp1q(tN,exp(.5*Lev),t1);
end
Nv = max(0,Nv); %in case negative values, set to 0

%Nm is the total modified loudness, Nm(t). This is LPFed by a 3-point
%triangular filter, n=26 times (I do 13 fwd and 13 backwd for zero-phase),
%which results in ~Gaussian smoothing with sig~=100 ms.
gv = ones([nczs 1],'double'); %weights
gv(czs<3) = 0; gv(czs>19) = -1;
Nm = Nv*gv;
b = ones([3 1],'double')./3; n = 13;
for nn = 1:n, Nm = filtfilt(b,1,Nm); end



%REFERENCES:
%
%Bismarck Gv (1973). Vorschlag f�r ein einfaches Verfahren zur
%Klassifikation station�rer Sprachschalle. Acustica 28(3): 186-7.
%
%Zwicker E, Terhardt E, Paulus E (1979). Automatic speech recognition using
%psychoacoustic models. J Acoust Soc Am 65(2): 487-98.
%
%Ruske G, Schotola T (1978). An approach to speech recognition using
%syllabic decision units. ICASSP: IEEE. 3: 722-5.
%
%Ruske G, Schotola T (1981). The efficiency of demisyllable segmentation in
%the recognition of spoken words. ICASSP: IEEE. 6: 971-4
%
%Schotola T (1984). On the use of demisyllables in automatic word
%recognition. Speech Commun 3(1): 63-87.

end


%% landmark detection in envelope
function [allTS, varNames] = find_landmarks(TS,  onsOff, cleanup_flag, mpd, minRatio)

%% find discrete events in envelope

% make sure envelope is row vector
if ~isrow(TS)
    TS = TS';
end

% Bring max envelope to 0
TS = TS-max(TS)-0.01;

envFloor = ceil(min(TS)+1);

% first temporal derivative of TS
diff_loudness = [0 diff(TS)];

%% discrete loudness

% min
% [lmin, minloc] = findpeaks(-TS,'MinPeakDistance',mpd);
[lmin, minloc] = findpeaks(-TS,'MinPeakProminence',6);
minEnv = zeros(size(TS));
minEnv(minloc)=-lmin;

% Also add the time that envelope crosses the floor threshold
% (adding this because min peaks are found at beginning, not end, of 
% intertrial intervals for speech
if envFloor<80
    idxOnset = 429+find(TS(430:505)>envFloor,1,'first');
    minEnv(idxOnset) = TS(idxOnset);
end


% max
% [lmin, minloc] = findpeaks(TS,'MinPeakDistance',mpd);
[lmin, minloc] = findpeaks(TS,'MinPeakProminence',6);
peakEnv = zeros(size(TS));
peakEnv(minloc)=lmin;


% Check for cases where two mins without a max in between them 
Mins = find(minEnv);
Maxs = find(peakEnv);

for im = 2:numel(Mins)
    
    iMax = find(Maxs>Mins(im-1) & Maxs<Mins(im));
    
    mE = envFloor-1;
    if isempty(iMax)
        [mE,imE] = max(TS(Mins(im-1):Mins(im)));
    end
    
    % Include only if envelope exceeds a low threshold
    if mE>envFloor+2
        peakEnv(imE+Mins(im-1)-1) = mE;
    end
end


%% discrete delta loudness

% min
negloud = diff_loudness; negloud(negloud>0) = 0;
[lmin, minloc] = findpeaks(-negloud,'MinPeakDistance',mpd);
% [lmin, minloc] = findpeaks(-negloud,'MinPeakProminence',6);
decrRate = zeros(size(TS));
decrRate(minloc)=-lmin;

% max
posloud = diff_loudness; posloud(posloud<0) = 0;
[lmin, minloc] = findpeaks(posloud,'MinPeakDistance',mpd);
% [lmin, minloc] = findpeaks(posloud,'MinPeakProminence',6);
peakRate = zeros(size(TS));
peakRate(minloc)=lmin;

clear negloud posloud ;


%% ------- KP clean up 
% for each minEnv to maxEnv, keep the highest peakRate event

Mins = find(minEnv);
Maxs = find(peakEnv);
iprs = find(peakRate);

NewPeakRate = zeros(size(peakRate));

for im = 1:numel(Mins)
    
    % Find beginning and end of this level ramp
    t1 = Mins(im);
    t2 = Maxs(find(Maxs>t1,1,'first'));
    if isempty(t2)
        continue
    end
    
    % Find maximum peak of derivative
    thesePRs = iprs(iprs>t1 & iprs<t2);
    [maxRise,iRise] = max(peakRate(thesePRs));
    
    % Set all other peakRate events to keep
    NewPeakRate(thesePRs(iRise)) = maxRise;
    
end

% % % Plot to check
% % figure;
% % plot(TS,'k','LineWidth',2)
% % hold on
% % plot(find(peakRate),TS(1,peakRate~=0),'.r','MarkerSize',10)
% % plot(find(peakEnv),TS(1,peakEnv~=0),'.g','MarkerSize',20)
% % plot(find(minEnv),TS(1,minEnv~=0),'.c','MarkerSize',20)
% % 
% % plot(find(NewPeakRate),TS(NewPeakRate~=0),'*b','MarkerSize',20)


%% Finish storing output information

allTS = [TS; ...
    diff_loudness;...
    minEnv;...
    peakEnv;...
    decrRate;...
    NewPeakRate];

varNames = {'Loudness', 'dtLoudness', 'minenv', 'peakEnv', 'minRate', 'peakRate'};

end
