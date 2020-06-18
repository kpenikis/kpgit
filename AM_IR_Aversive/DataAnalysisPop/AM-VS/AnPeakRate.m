function AnPeakRate
%
%  AnPeakRate
%    Yulia analysis, compare AM and speech
%
%  KP, 2019-12
%

% close all
global fn spkshift smth_win exclOnset AM_durs VS_durs nTrGrp

%!!!!!!!!!!!!!!!!!!!!
smth_win    = 10;
exclOnset   = 0;
%!!!!!!!!!!!!!!!!!!!!
nTrGrp = 4;

histedges = logspace(-4,-1,100);


%% Units setup

fn = set_paths_directories('','',1);

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
        %         keyboard
    end
end



%% Data setup

AMrates  = [2 4 8 16 32];
AM_durs  = [1500   1000   1000   1000   1000   1000   1937   1937   1000];
VS_durs  = [1670   1406   2412   5928   2554   2556];


%% Figure settings

savedir = fullfile(fn.figs,'TRF_AMVS');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallrect    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];

dotSize = 40;

% Set colors
subjcolors = [0.5 0.2 0.7];
Subjects   = unique({UData_AM.Subject});


%%

% Begin by loading just one session (Mar30-AM and Mar30-VS)
allTS_AM=[];
allTS_Sp=[];

for iiUn = 1:size(Uindices_AMVS,1)
    
    iUnAM = Uindices_AMVS(iiUn,1);
    iUnVS = Uindices_AMVS(iiUn,2);
    
    if isnan(iUnAM)
        continue
    end
    
    % Get PSTH and STIM data
    
    [PSTH_AMc,STIM_AMc] = getdata4TRF(UData_AM,iUnAM);
    [PSTH_Spc,STIM_Spc] = getdata4TRF(UData_VS,iUnVS);
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Concatenate AM // Sp data
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    psth_AM = []; stim_AM = [];    
    for st = find(~cellfun(@isempty,PSTH_AMc))'
        if size(STIM_AMc{st},1)<10
            keyboard
        end
%         if st==5
%             continue
%         end
        
        % Average trials and concatenate stim
        psth_AM  = [psth_AM mean(PSTH_AMc{st}(:,:),1) ];
        stim_AM  = [stim_AM mean(STIM_AMc{st}(:,:),1) ];
    end
    
    psth_Sp = []; stim_Sp = [];
    for st = find(~cellfun(@isempty,PSTH_Spc))'
        if size(STIM_Spc{st},1)<10
            continue
        end
        
        psth_Sp  = [psth_Sp mean(PSTH_Spc{st}(:,:),1) ];
        stim_Sp  = [stim_Sp mean(STIM_Spc{st}(:,:),1) ];
        
    end
    
    % Normalize both stim signals to be in range [0 1]
%     stim_AM = stim_AM/max(stim_AM);
%     stim_Sp = stim_Sp/max(stim_Sp);
    
    
    
    % Correct artifact at beginning
    [minV,iminV] = min(stim_AM(1:30));
    stim_AM(1:iminV) = minV;
    
    % Find sound events
%     allTS = [TS; ...
%              diff_loudness;...
%              minEnv;...
%              peakEnv;...
%              minRate;...
%              peakRate];
    uTS_AM = find_peakRate(stim_AM, 1000, [1 length(stim_AM)]./1000, 'broadband');
    uTS_Sp = find_peakRate(stim_Sp, 1000, [1 length(stim_Sp)]./1000, 'broadband');
    
    allTS_AM = [allTS_AM uTS_AM];
    allTS_Sp = [allTS_Sp uTS_Sp];
    
    % Save data from this unit for plotting
    
    % actually, maybe just stim and psth
    % then average across all units
    % careful that stimuli line up correctly: save AC and DB and both
    % separately, or insert nans and ignore
    % also: different spl & lp values
    
    % histcounts for peakRate and peakEnv
    
    % psth for later averaging and calculating 
    
    
end %iiUn


% allTS = [TS; ...
%          diff_loudness;...
%          minEnv;...
%          peakEnv;...
%          minRate;...
%          peakRate];

% a) Plot distribution of peakRate events in AM and speech
% b) Evoked response from all cells averaged (color by peakRate magnitude)
% c) Evoked response for each cell individually (mark t_peak)
% d) Distribution of t_peaks across cells

hf=figure;

% a) histograms of peakRate
subplot(4,2,1:2)
hold on
histogram(allTS_Sp(6,:),histedges,'FaceColor','m')
histogram(allTS_AM(6,:),histedges,'FaceColor','k')
set(gca,'xscale','log','xlim',10.^[-4 -1])
xlabel('peakRate magnitude')
ylabel('Count')

% Set and save Boundaries for peakRate


% b)




end