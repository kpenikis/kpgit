function LatDur_CTTS
% 
% PopResp_CTTS
% 

close all


whichClass   = 'Full';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELLS
whichCells   = 'SpecRest'; %'Q_pkFR'; %'dpRank_RS'; % pkFR_RS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME
Dur          = 500;
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
% flagRS = find(UnitInfo(theseUns,:).TroughPeak>0.43);
flagRS = UnitInfo(theseUns,:).TroughPeak>0.43;
flagNS = UnitInfo(theseUns,:).TroughPeak<=0.43;


% Flag non-significant SU dps
UnSig = bootstrap4significance(CReach);


FRdiffThresh = 5;
% pdms = 500;
DurMax = 500;
upd_lat = makedist('uniform','Lower',1,'Upper',500);
upd_dur = makedist('uniform','Lower',1,'Upper',DurMax);


if strcmp(whichCells,'SpecRest')
    
    switch whichStim
        case 'AC'
            iC_Rest = find(CReach.dprime<1 & flagRS & UnSig);
            iC_Spec = find(CReach.dprime>1 & flagRS);
            iC_NS   = find(flagNS);
        case 'Speech'
            iC_Rest = find(CReach.dprime<1.4 & flagRS & UnSig);
            iC_Spec = find(CReach.dprime>1.4 & flagRS);
            iC_NS   = find(flagNS);
    end
    
    % [ Sp/Rest lat dur d' ]
    Result = nan(size(CReach,1),5,8);
    ResultStats = struct;
    
    for ist = 1:8
        
        hf(ist)=figure;
        set(gcf,'Position',fullscreen)
        hold on
        
        %------------------------------------------------------------------
        %%----Rest (RS non-specialist cells)
        %------------------------------------------------------------------
        
        for iu = 1:numel(iC_Rest)
            
            pk_time   = nan;
            t_onset   = nan;
            t_offset  = nan;
            t_mid     = nan;
            t_up      = nan;
            
            % Get latency and duration
            Data = permute(mean(CTTS(iC_Rest(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
            
            % Find half max FR and peak events above it
            peakThresh = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /4) + min(Data(:,AnWin));
            if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
                continue
            end
            [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',peakThresh);
            
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
                if Data(:,ims)<peakThresh
                    t_onset = ims - AnWin(1) - 1 ;
                    break
                end
            end
            
            % Find time just after falling below halfMax
            %         for ims = [pk_time:pdms 1:pk_time]
            for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2)
                if Data(:,ims)<peakThresh
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
            Result(iC_Rest(iu),1,ist)  = 1;
            
            %change in 4 places
            thisLat = pk_time; %t_onset; %t_mid;
            
            Result(iC_Rest(iu),2,ist)  = thisLat;
            Result(iC_Rest(iu),3,ist)  = t_up;
            Result(iC_Rest(iu),4,ist)  = ceil(pk_height);
            Result(iC_Rest(iu),5,ist)  = CReach.dprime(iC_Rest(iu));
            
            % Lat vs dur (size peak height)
            subplot(4,4,[9 10 13 14]);
            hold on
            scatter(thisLat,t_up,5*ceil(pk_height),'o','MarkerFaceColor','k','MarkerEdgeColor','none')
            
        end %iu
        
        subplot(4,4,5:6)
        hold on
        histogram(Result(iC_Rest,2,ist),1:10:500,'FaceColor','k','EdgeColor','none')
        
        subplot(4,4,[11 15])
        hold on
        histogram(Result(iC_Rest,3,ist),1:10:500,'FaceColor','k','EdgeColor','none')
        
        
        %------------------------------------------------------------------
        %%%----Specialists
        %------------------------------------------------------------------
        
        for iu = 1:numel(iC_Spec)
            
            pk_time   = nan;
            t_onset   = nan;
            t_offset  = nan;
            t_mid     = nan;
            t_up      = nan;
            
            % Get latency and duration
            Data = permute(mean(CTTS(iC_Spec(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
            
            % Find half max FR and peak events above it
            peakThresh = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /4) + min(Data(:,AnWin));
            if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
                continue
            end
            [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',peakThresh);
            
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
                if Data(:,ims)<peakThresh
                    t_onset = ims - AnWin(1) - 1 ;
                    break
                end
            end
            
            % Find time just after falling below halfMax
            %         for ims = [pk_time:pdms 1:pk_time]
            for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2)
                if Data(:,ims)<peakThresh
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
            Result(iC_Spec(iu),1,ist)  = 1;
            
            thisLat = pk_time; %t_onset; %t_mid;
            
            Result(iC_Spec(iu),2,ist)  = thisLat;
            Result(iC_Spec(iu),3,ist)  = t_up;
            Result(iC_Spec(iu),4,ist)  = ceil(pk_height);
            Result(iC_Spec(iu),5,ist)  = CReach.dprime(iC_Spec(iu));
            
            % Lat vs dur (size peak height)
            subplot(4,4,[9 10 13 14]);
            hold on
            scatter(thisLat,t_up,5*ceil(pk_height),'o','MarkerFaceColor','k','MarkerEdgeColor','none')

        end %iu
        
        subplot(4,4,5:6)
        hold on
        histogram(Result(iC_Spec,2,ist),1:10:500,'FaceColor','k','EdgeColor','none')
        
        subplot(4,4,[11 15])
        hold on
        histogram(Result(iC_Spec,3,ist),1:10:500,'FaceColor','k','EdgeColor','none')
        
        
        
        %------------------------------------------------------------------
        %%----NS cells
        %------------------------------------------------------------------
        
        for iu = 1:numel(iC_NS)
            
            pk_time   = nan;
            t_onset   = nan;
            t_offset  = nan;
            t_mid     = nan;
            t_up      = nan;
            
            % Get latency and duration
            Data = permute(mean(CTTS(iC_NS(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
            
            % Find half max FR and peak events above it
            peakThresh = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /4) + min(Data(:,AnWin));
            if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
                continue
            end
            [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',peakThresh);
            
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
                if Data(:,ims)<peakThresh
                    t_onset = ims - AnWin(1) - 1 ;
                    break
                end
            end
            
            % Find time just after falling below halfMax
            %         for ims = [pk_time:pdms 1:pk_time]
            for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2)
                if Data(:,ims)<peakThresh
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
            Result(iC_NS(iu),1,ist)  = 1;
            
            %change in 4 places
            thisLat = pk_time; %t_onset; %t_mid;
            
            Result(iC_NS(iu),2,ist)  = thisLat;
            Result(iC_NS(iu),3,ist)  = t_up;
            Result(iC_NS(iu),4,ist)  = ceil(pk_height);
            Result(iC_NS(iu),5,ist)  = CReach.dprime(iC_NS(iu));
            
            % Lat vs dur (size peak height)
            subplot(4,4,[9 10 13 14]);
            hold on
            scatter(thisLat,t_up,5*ceil(pk_height),'o','MarkerFaceColor','g','MarkerEdgeColor','none')
            
        end %iu
        
        subplot(4,4,5:6)
        hold on
        histogram(Result(iC_NS,2,ist),1:10:500,'FaceColor','g','EdgeColor','none')
        
        subplot(4,4,[11 15])
        hold on
        histogram(Result(iC_NS,3,ist),1:10:500,'FaceColor','g','EdgeColor','none')
        
        
        
        %------------------------------------------------------------------
        % Finish plots
        %------------------------------------------------------------------
        
        subplot(4,4,[9 10 13 14]);
        hold on
        xlim([0 500])
        ylim([0 DurMax])
        set(gca,'Color','none')
        xlabel('time exceeds 1/4 max')
        ylabel('time above 1/4 max')
        
        subplot(4,4,5:6) 
        hold on
        xlim([0 500])
        set(gca,'Color','none')
        
        subplot(4,4,[11 15])
        hold on
        xlim([0 DurMax])
        set(gca,'Xdir','reverse','Color','none')
        camroll(270)
        
        
        [h_un_lat,p_un_lat] = kstest(Result([iC_Rest; iC_Spec],2,ist),'cdf',upd_lat);
        [h_un_dur,p_un_dur] = kstest(Result([iC_Rest; iC_Spec],3,ist),'cdf',upd_dur);
        
        ResultStats.StimLat_Uniform(ist) = p_un_lat;
        ResultStats.StimDur_Uniform(ist) = p_un_dur;
        
        p_un_lat
        if h_un_lat==0
            fprintf('latency distr uniform for stim %i\n',ist)
            
        end
%         if h_un_dur==0
%             fprintf('duration distr uniform for stim %i\n',ist)
%             p_un_dur
%         end
        
        suptitle(num2str(ist))
        
        
        savename = sprintf('%s_%s_stim%i_%s',whichStim,whichCells,ist,'pk_time_wNS');
%         print_eps_kp(hf(ist),fullfile(figsavedir,savename))
        
    end %ist
    
    keyboard
    % Finish up (stats)
    
    Lats_Rest = permute(Result(iC_Rest,2,:),[1 3 2]);
    Lats_Spec = permute(Result(iC_Spec,2,:),[1 3 2]);
    
    % Bootstrap to test for difference in distributions
    nbs = 10000;
    pvals   = nan(nbs,1);
    medians = nan(nbs,1);
    for ibs=1:nbs
        irnd_Spec = ceil(rand(size(Lats_Spec,1),1)*size(Lats_Spec,1));
        irnd_Rest = ceil(rand(size(Lats_Rest,1),1)*size(Lats_Rest,1));
        data_Spec = Lats_Spec(irnd_Spec,:);
        data_Rest = Lats_Rest(irnd_Rest(1:length(irnd_Spec)),:);
        [~,pvals(ibs)]=kstest2(data_Spec(:),data_Rest(:));
        
        medians(ibs) = median(data_Spec(:),'omitnan') - median(data_Rest(:),'omitnan');
    end
    figure;
    histogram(medians)
    sum(pvals<0.01)
    
    [h,p]=kstest2(Lats_Rest(:),Lats_Spec(:));
    fprintf('\nLatency distributions for Specialists and Rest:\n  2-sample kstest p=%0.1e\n',p)
    
    ResultStats.Lat_SpecRest_kstest2 = p;
    
    Durs_Rest = permute(Result(iC_Rest,3,:),[1 3 2]);
    Durs_Spec = permute(Result(iC_Spec,3,:),[1 3 2]);
    [h,p]=kstest2(Durs_Rest(:),Durs_Spec(:));
    fprintf('\nDuration distributions for Specialists and Rest:\n  2-sample kstest p=%0.1e\n\n',p)
    
    ResultStats.Dur_SpecRest_kstest2 = p;
    
    savename = sprintf('%s_%s_%s',whichStim,whichCells,'pk_time');
    save(fullfile(figsavedir,savename),'ResultStats','Result','-v7.3')
end


keyboard


if strcmp(whichCells,'Q_pkFR')
    
    % [ Q lat dur d' ]
    Result = nan(sum(flagRS),4,8);
    
    for ist = 1:8
        
        hf(ist)=figure;
        set(gcf,'Position',fullscreen)
        
        
        ii=0;
        % For each Quantile
        for iq = 1:size(CR_Qpfr_Sp,1)
            
            % For each SU in this peakFR Quantile
            iC_Spec = CR_Qpfr_Sp(iq,:).CRids{:}; %only RS cells
            
            for iu = 1:numel(iC_Spec)
                
                pk_time   = nan;
                t_onset   = nan;
                t_offset  = nan;
                t_mid     = nan;
                t_up      = nan;
                
                % Get latency and duration
                Data = permute(mean(CTTS(iC_Spec(iu),:,:,ist),3,'omitnan'),[1 2 4 3]).*1000;
                
                % Find half max FR and peak events above it
                peakThresh = (( max(Data(:,AnWin)) - min(Data(:,AnWin)) ) /2) + min(Data(:,AnWin));
                if max(Data(:,AnWin))<FRdiffThresh/5 || (max(Data(:,AnWin))-min(Data(:,AnWin))) < FRdiffThresh
                    continue
                end
                [PKS,LOCS] = findpeaks(Data(:,AnWin),'MinPeakHeight',peakThresh);
                
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
                    if Data(:,ims)<peakThresh
                        t_onset = ims - AnWin(1) - 1 ;
                        break
                    end
                end
                
                % Find time just after falling below halfMax
                %         for ims = [pk_time:pdms 1:pk_time]
                for ims = ( AnWin(1) - 1 + pk_time ) : size(Data,2)
                    if Data(:,ims)<peakThresh
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
                Result(ii,4,ist)  = CReach.dprime(iC_Spec(iu));
                
                
                % PLOTS
                
                % Lat vs dur (size d')
                subplot(3,size(CR_Qpfr_Sp,1),iq);
                hold on
                plot(thisLat,t_up,'.','Color','b','MarkerSize',20*(0.5+CReach.dprime(iC_Spec(iu))) )
                xlim([0 500])
                ylim([0 500])
                
                % Lat vs d'
                subplot(3,size(CR_Qpfr_Sp,1),iq+5);
                hold on
                plot(thisLat,CReach.dprime(iC_Spec(iu)),'.','Color','k','MarkerSize',10)
                xlim([0 500])
                ylim([-0.5 3])
                
                % Dur vs d'
                subplot(3,size(CR_Qpfr_Sp,1),iq+10);
                hold on
                plot(t_up,CReach.dprime(iC_Spec(iu)),'.','Color','k','MarkerSize',10)
                xlim([0 500])
                ylim([-0.5 3])
                
            end
        end
        
        suptitle(['st ' num2str(ist)])
        
        savename = sprintf('%s_stim%i_%s_%s',whichStim,ist,whichCells,whichClass);
        print_eps_kp(gcf,fullfile(figsavedir,savename))
    end
    
end





end





