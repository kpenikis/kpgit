function PowerSpectra_TRFcode
%
%  PowerSpectra_TRFcode
%
%  KP, 2020-05
%



close all
global fn spkshift smth_win exclOnset AM_durs VS_durs nTrGrp

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
smth_win    = 5;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
exclOnset   = 0; 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nTrGrp      = 10;


%% Units setup

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UData_AM = q.UnitData;
UInfo_AM = q.UnitInfo;
clear q
%-------
spkshift = 0; %mean([UData_AM([UData_AM.IntTime_spk]>0).IntTime_spk]);
%-------

q = load(fullfile(fn.processed,'UnitsVS'));
UData_VS = q.UnitData;
UInfo_VS = q.UnitInfo;
clear q

% Load AM-VS LUT
load(fullfile(fn.processed,'UnMatchedLUT'));

Uindices_AMVS = nan(size(MatchedUnits,1),2);
for iu = 1:size(MatchedUnits,1)
    iU_AM = []; iU_VS = []; 
    iU_AM = find( strcmp(MatchedUnits.AMsessions{iu},{UData_AM.Session}) & MatchedUnits.AMclu(iu)==[UData_AM.Clu] );
    iU_VS = find( strcmp(MatchedUnits.VSsessions{iu},{UData_VS.Session}) & MatchedUnits.VSclu(iu)==[UData_VS.Clu] );
    if numel(iU_AM)==1 && numel(iU_VS)==1
        Uindices_AMVS(iu,:) = [iU_AM iU_VS];
    else
        keyboard
    end
end


%% Data setup

AMrates   = [2 4 8 16 32];
AM_durs = [1500   1000   1000   1000   1000   1000   1937   1937   1000];
VS_durs = [1670   1406   2412   5928   2554   2556];



%% Figure settings

savedir = fullfile(fn.figs,'StimModPower');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',14)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallrect    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];


%% Collect data

STIM_AM_ALL = cell(7,1);
STIM_VS_ALL = cell(6,1);

for iiUn = 1:size(Uindices_AMVS,1)
    
    iUnAM = Uindices_AMVS(iiUn,1);
    iUnVS = Uindices_AMVS(iiUn,2);
    
    
    % Get PSTH and Stim data
    [~,STIM_AM] = getdata4TRF(UData_AM,iUnAM);
    [~,STIM_VS] = getdata4TRF(UData_VS,iUnVS);
    
    % Save averages for fft later
    for ist = 1:numel(STIM_AM)
        STIM_AM_ALL{ist} = [STIM_AM_ALL{ist}; mean(STIM_AM{ist},1)];
    end
    for ist = 1:numel(STIM_VS)
        STIM_VS_ALL{ist} = [STIM_VS_ALL{ist}; mean(STIM_VS{ist},1)];
    end
    
end %iUn



%% Plot MPS

fs = 1000;


% Modulation Power Spectrum of all stimuli
CatEnv_Pdc   = [];
for ist =  1:5 %numel(STIM_AM_ALL) 
    CatEnv_Pdc  = [CatEnv_Pdc;  mean(STIM_AM_ALL{ist},1)'];
end
CatEnv_Irr   = [];
for ist =  [6 7]
    CatEnv_Irr  = [CatEnv_Irr;  mean(STIM_AM_ALL{ist},1)'];
end
CatEnv_VS   = [];
for ist =  1:numel(STIM_VS_ALL) 
    CatEnv_VS  = [CatEnv_VS; zeros(200,1); mean(STIM_VS_ALL{ist},1)'];
end

% FFT method #1
% [ssf,fws]     = plotMPS(CatEnv,fs);
% [ssf_R,fws_R] = plotMPS(CatPSTH,fs);

% FFT method #2

% AM - Pdc
clear f y_S P1_S P2_S

L    = length(CatEnv_Pdc);
% f    = logspace(-3,2.699,floor(L/2+1))';
f_Pd    = fs*(0:(L/2))'/L;
y_S  = fft(CatEnv_Pdc);
P2_S = abs(y_S/L);
P1_S = P2_S( 1 : floor(L/2+1) );
P1_S = P1_S.*sqrt(f_Pd);
P1_S(2:end-1) = 2*P1_S(2:end-1);

MPS_Pd = P1_S;
MPS_Pd = conv(MPS_Pd,gausswin(3),'same');
f_Pd   = f_Pd(unique(round(logspace(0,log10(L/2),100))));
MPS_Pd = MPS_Pd(unique(round(logspace(0,log10(L/2),100))));
MPS_Pd = MPS_Pd./max(MPS_Pd);



% AM - Irr
clear f y_S P1_S P2_S

L    = length(CatEnv_Irr);
% f    = logspace(-3,2.699,floor(L/2+1))';
f_AM    = fs*(0:(L/2))'/L;
y_S  = fft(CatEnv_Irr);
P2_S = abs(y_S/L);
P1_S = P2_S( 1 : floor(L/2+1) );
P1_S = P1_S.*sqrt(f_AM);
P1_S(2:end-1) = 2*P1_S(2:end-1);

MPS_AM = P1_S;
MPS_AM = conv(MPS_AM,gausswin(3),'same');
f_AM   = f_AM(unique(round(logspace(0,log10(L/2),100))));
MPS_AM = MPS_AM(unique(round(logspace(0,log10(L/2),100))));
MPS_AM = MPS_AM./max(MPS_AM);


% VS
clear f y_S P1_S P2_S

L    = length(CatEnv_VS);
f_VS    = fs*(0:(L/2))'/L;
y_S  = fft(CatEnv_VS);
P2_S = abs(y_S/L);
P1_S = P2_S( 1 : floor(L/2+1) );
P1_S = P1_S.*sqrt(f_VS);
P1_S(2:end-1) = 2*P1_S(2:end-1);

MPS_VS = P1_S;
MPS_VS = conv(MPS_VS,gausswin(3),'same');
f_VS   = f_VS(unique(round(logspace(0,log10(L/2),100))));
MPS_VS = MPS_VS(unique(round(logspace(0,log10(L/2),100))));
MPS_VS = MPS_VS./max(MPS_VS);

% Plot SIGNAL
% figure; 
% plot(CatEnv_AM./max(CatEnv_AM),'LineWidth',2)
% xlim([0 length(CatEnv_AM)])
% xlabel('Time (ms)')

% print_eps_kp(gcf,fullfile(savedir,'Data',['CatStimPSTH_' whichStim '_10msG']))


% Plot POWER SPECTRUM

figure; 
xlim([0.25 2^7])
set(gca,'xscale','log')
% plot(ssf,fws./max(fws(ssf>0.25&ssf<64)),'LineWidth',2)
hold on
plot(f_Pd,MPS_Pd,'LineWidth',2)
plot(f_AM,MPS_AM,'LineWidth',2)
plot(f_VS,MPS_VS,'LineWidth',2)
set(gca,'xtick',[0.5 1 2 4 8 16 32 2^7],'Color','none')
xlabel('Frequency')

legend({'Pdc' 'Irr' 'Speech'})

print_eps_kp(gcf,fullfile(savedir,'MPS_SinSpeech_2')) %'_10msG'


%% Speech: each stim MPS separaetely, then average


MPS_VS_all = [];

for ist =  1:numel(STIM_VS_ALL)
    CatEnv_VS  = mean(STIM_VS_ALL{ist},1)';
    figure;
    dwt(CatEnv_VS,1000)
    clear f y_S P1_S P2_S
    
    L    = length(CatEnv_VS);
    f_VS    = fs*(0:(L/2))'/L;
    
    y_S  = fft(CatEnv_VS);
    

% % n = length(CatEnv_VS);          % number of samples
% % f = (0:n-1)*(fs/n);     % frequency range
% % power = abs(y).^2/n;    % power of the DFT
% % 
% % figure;
% % plot(f,power)
% % xlabel('Frequency')
% % ylabel('Power')
% % xlim([0 20])

    P2_S = abs(y_S/L);
    P1_S = P2_S( 1 : floor(L/2+1) );
    P1_S = P1_S.*sqrt(f_VS);
    P1_S(2:end-1) = 2*P1_S(2:end-1);
    
    MPS_VS = P1_S;
    MPS_VS = conv(MPS_VS,gausswin(3),'same');
    f_VS   = f_VS(unique(round(logspace(0,log10(L/2),100))));
    MPS_VS = MPS_VS(unique(round(logspace(0,log10(L/2),100))));
%     MPS_VS = MPS_VS./max(MPS_VS);
    
    MPS_VS_all = [MPS_VS_all MPS_VS];
end



figure; 
pwelch(STIM_VS_ALL{ist}')



end
    