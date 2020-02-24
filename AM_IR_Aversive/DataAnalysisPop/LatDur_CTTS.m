function LatDur_CTTS
% 
% PopResp_CTTS
% 

close all


whichClass   = 'Full';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'Q_pkFR'; %'dpRank_RS'; % pkFR_RS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 300;
WinBeg       = 501 * ones(size(Dur));
WinEnds      = WinBeg+Dur-1;
AnWin        = WinBeg:WinEnds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIM
whichStim    = 'Speech';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolution
tau          = 5;
lambda       = 1/tau;
% winlen       = 500;
convwin      = exp(-lambda*(1:500));
convwin      = convwin./sum(convwin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        rootdir = fullfile(fn.figs,'ClassAM');
        rawdata = 'CTTS_AM';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'Units'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassSpeech');
        rawdata = 'CTTS_Speech_nonSim';
        
        % Load Unit data files
        q = load(fullfile(fn.processed,'UnitsVS'));
        UnitData = q.UnitData;
        UnitInfo = q.UnitInfo;
        clear q
        
end

% Load spikes data (created in gatherCellTimeTrialStim, used to be cumulativeSpikeCount)
q=load(fullfile(rootdir,'RawData',rawdata)); %Cell_Time_Trial_Stim
Cell_Time_Trial_Stim = q.Cell_Time_Trial_Stim;
Env_Time_Trial_Stim  = q.Env_Time_Trial_Stim;

% Load SU classification results
q = load(fullfile(rootdir,whichStim,whichClass,'each','CR_each.mat'));
CReach = q.CR;
clear q

% Load Quantile classification results
q=load(fullfile(rootdir,whichStim,whichClass,'Q_pkFR',['CR_v' whichClass '_Q_pkFR.mat']));
CR_Qpfr_Sp = q.CR;
clear q

% Load encoding time estimation
load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassSpeech/avgPropUp.mat')


%%
% Figure settings
set(groot,'DefaultTextInterpreter','none')
set(groot,'DefaultAxesFontSize',18)
set(groot,'defaultAxesTickDir', 'out');
set(groot,'defaultAxesTickDirMode', 'manual');

scrsz = get(0,'ScreenSize');     %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];
tallhalf    = [1 scrsz(4) scrsz(3)/2 scrsz(4)];
widehalf    = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

% Set figsavedir
figsavedir = fullfile(fn.figs,'PopResp','LatDur');
if ~exist(figsavedir,'dir')
    mkdir(figsavedir)
end


%% Prepare to parse data

nTrialMat = nan(size(Cell_Time_Trial_Stim,1),size(Cell_Time_Trial_Stim,4));
for ist = 1:size(Cell_Time_Trial_Stim,4)
    CT  = permute(sum(Cell_Time_Trial_Stim(:,:,:,ist),2),[1 3 2]);
    nTrialMat(:,ist) = sum(~isnan(CT),2);
end

switch whichStim
    case 'AC'
        theseStim  = 1:8;
    case 'DB'
        theseStim  = [1:6 9:10];
    case 'Speech'
        theseStim  = [4 3 2 1 5 6 7 8]; %1:size(Cell_Time_Trial_Stim,4);
end
% theseStim = theseStim(1:6);

% CellTypes
flagRS = find(UnitInfo.TroughPeak>0.43);
flagNS = find(UnitInfo.TroughPeak<=0.43);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
minTrs = 12;

% Get UnitData indices for matching to CR each
[CTTS,theseUns,~,~,~] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    'each', nTrialMat, UnitData,theseStim, flagRS, flagNS, minTrs, convwin, 1:size(Cell_Time_Trial_Stim,2), 'exp' );

if size(CTTS,1)~=size(CReach,1)
    keyboard
end

% New RS/NS labels
flagRS = find(UnitInfo(theseUns,:).TroughPeak>0.43);


FRdiffThresh = 5;
% pdms = 500;



% [ Q lat dur d' ]
Result = nan(numel(flagRS),4,8);

for ist = 1:8
    
    hf(ist)=figure;
    set(gcf,'Position',fullscreen)
    
    
    ii=0;
    % For each Quantile
    for iq = 1:size(CR_Qpfr_Sp,1)
        
        % For each SU in this peakFR Quantile
        iCRe = CR_Qpfr_Sp(iq,:).CRids{:}; %only RS cells
        
        for iu = 1:numel(iCRe)
            
            pk_time   = nan;
            t_onset   = nan;
            t_offset  = nan;
            t_mid     = nan;
            t_up      = nan;
            
            % Get latency and duration
            Data = permute(mean(CTTS(iCRe(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
            
            % Find half max FR and peak events above it
            halfMax = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /2) + min(Data(:,AnWin));
            if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
                continue
            end
            [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',halfMax);
            
            if isempty(PKS)
                continue
            end
            
            % ASSUME just one peak
            [~,ipk] = max(PKS);
            pk_height = PKS(ipk);
            pk_time   = LOCS(ipk);
            
            
            % Find time just before surpassing halfMax
            %         for ims = [pk_time:-1:1 pdms:-1:pk_time]
            for ims = ( AnWin(1) - 1 +  pk_time ) : -1 : 1
                if Data(:,ims)<halfMax
                    t_onset = ims - AnWin(1) - 1 ;
                    break
                end
            end
            
            % Find time just after falling below halfMax
            %         for ims = [pk_time:pdms 1:pk_time]
            for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2) 
                if Data(:,ims)<halfMax
                    t_offset = ims - AnWin(1) - 1 ;
                    break
                end
            end
            
            % Calculate peak duration
            if t_offset>=t_onset
                t_up = t_offset-t_onset;
%             else
%                 t_up = t_offset+pdms-t_onset;
            end
            
            % Time of midpoint of peak
            t_mid = t_onset+t_up/2;
%             if t_mid>pdms
%                 t_mid = t_mid-pdms;
%             end
            
            
            % Save data for statistics
            ii = ii+1;
            Result(ii,1,ist)  = iq;
            thisLat = t_mid; %t_onset; %pk_time; %
            Result(ii,2,ist)  = thisLat;
            Result(ii,3,ist)  = t_up;
            Result(ii,4,ist)  = CReach.dprime(iCRe(iu));
            
            
            % PLOTS
            
            % Lat vs dur (size d')
            subplot(3,size(CR_Qpfr_Sp,1),iq);
            hold on
            plot(thisLat,t_up,'.','Color','b','MarkerSize',20*(0.5+CReach.dprime(iCRe(iu))) )
            xlim([0 500])
            ylim([0 500])
            
            % Lat vs d'
            subplot(3,size(CR_Qpfr_Sp,1),iq+5);
            hold on
            plot(thisLat,CReach.dprime(iCRe(iu)),'.','Color','k','MarkerSize',10)
            xlim([0 500])
            ylim([-0.5 3])
            
            % Dur vs d'
            subplot(3,size(CR_Qpfr_Sp,1),iq+10);
            hold on
            plot(t_up,CReach.dprime(iCRe(iu)),'.','Color','k','MarkerSize',10)
            xlim([0 500])
            ylim([-0.5 3])
            
        end
    end
    
    suptitle(['st ' num2str(ist)])
    
    savename = sprintf('%s_stim%i_%s_%s',whichStim,ist,whichCells,whichClass);
    print_eps_kp(gcf,fullfile(figsavedir,savename))
end







end





