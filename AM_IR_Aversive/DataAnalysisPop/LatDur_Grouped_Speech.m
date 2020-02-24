
% close all
useFR    =   'log'; 
sortStim = 3;

% Load Unit data files
fn = set_paths_directories('','',1);
q = load(fullfile(fn.processed,'UnitsVS'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q


% Load RepSpeech templates
k=load(fullfile(fn.stim,'SpeechStim','RepeatedSpeechTemplates'));
kfns = fieldnames(k);

Duration = max(structfun(@length,k));
Segs     = 1:numel(structfun(@length,k));
SegDurs  = structfun(@length,k);


savedir = fullfile(fn.figs,'PopulationSpeechSegments');
load(fullfile(savedir,'RepSpeech'))

if size(UnitInfo,1) ~= size(FR_vec,1)
    keyboard
end


%% Set up

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',18)

scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen  = [1 scrsz(4) scrsz(3) scrsz(4)];

colors = [   0   0   0;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
       
%% Latency  X  Duration

% Data settings
switch useFR
    case 'z'
        ThisData     = zFR_vec;
        Boundaries   = [-1 0 0.25 0.5 1 2];
        FRdiffThresh = 0.5;
    case 'log'
        ThisData     = FR_vec;
        Boundaries   = [-1 round(10.^[0.5 1 1.25 1.5]) ];
        FRdiffThresh = 4;
%         Boundaries  = [-1 round(10.^[0.5]) ];
end


%% Establish peak FR groups

%%%%  determine group id once, in classifier, 
%%%%  then apply to pop psths and this 
%%%%  *diff # cells
%%%%  OR just use same method (mean of max of each stim)

% CellTypes
flagRS = find(UnitInfo.TroughPeak>0.43);
flagNS = find(UnitInfo.TroughPeak<=0.43);

% RS cells
[pkRS,~]  = rankPeakFR(ThisData(flagRS,:,:,:));    % to group
[~,ipkRS] = max(ThisData(flagRS,:,sortStim),[],2); % to sort
Qbounds      = quantile(pkRS,[0 0.2 0.4 0.6 0.8 1]);

dataRS  = [];
isortRS = [];
for iq=1:5
    idx = find(pkRS>=Qbounds(iq) & pkRS<=Qbounds(iq+1));
    [lats,idx_sort] = sort(ipkRS(idx),'descend');
    isortRS  = [isortRS; idx(idx_sort)];
    dataRS   = [dataRS; pkRS(idx(idx_sort)) ipkRS(idx(idx_sort))];
%     ytickset = [ytickset size(lats,1)];
%     yticklab = [yticklab ['RS' num2str(iq)]];
end


% NS cells
[pkNS,~]  = rankPeakFR(ThisData(flagNS,:,:,:));    % to group
[~,ipkNS] = max(ThisData(flagNS,:,sortStim),[],2); % to sort
Qbounds      = quantile(pkNS,[0 1]);

dataNS  = [];
isortNS = [];
for iq=1
    idx = find(pkNS>=Qbounds(iq) & pkNS<=Qbounds(iq+1));
    [lats,idx_sort] = sort(ipkNS(idx),'descend');
    isortNS  = [isortNS; idx(idx_sort)];
    dataNS   = [dataNS; pkNS(idx(idx_sort)) ipkNS(idx(idx_sort))];
    ytickset = [ytickset size(lats,1)];
    yticklab = [yticklab 'NS'];
end

sortdata     = [dataRS; dataNS];
i_sorted     = [flagRS(isortRS); flagNS(isortNS)];
ytickset     = cumsum(ytickset);



%% Max Peak: time x duration


    
hfPk=figure; 
set(hfPk,'Position',fullscreen)

% show non sync units?
% plot cumulative nspikes behind? 

t_onset   = nan(size(ThisData,1),size(ThisData,3));
t_offset  = nan(size(ThisData,1),size(ThisData,3));
t_up      = nan(size(ThisData,1),size(ThisData,3));
t_mid     = nan(size(ThisData,1),size(ThisData,3));
pk_height = nan(size(ThisData,1),size(ThisData,3));
pk_time   = nan(size(ThisData,1),size(ThisData,3));

for iUn = 1:size(ThisData,1)
    
    for ist = 1:size(ThisData,3)
        
%         plot(ThisData(iUn,:,ist),'Color',colors(ist,:),'LineWidth',3)
%         hold on
        
        pdms = SegDurs(ist);
        
        halfMax = (( max(ThisData(iUn,:,ist)) - min(ThisData(iUn,:,ist)) ) /2) + min(ThisData(iUn,:,ist));
        if max(ThisData(iUn,:,ist))<FRdiffThresh/5 || (max(ThisData(iUn,:,ist))-min(ThisData(iUn,:,ist))) < FRdiffThresh
            continue
        end
        [PKS,LOCS] = findpeaks(ThisData(iUn,:,ist),'MinPeakHeight',halfMax);
        
        if isempty(PKS)
            continue
        end
        
        % ASSUME just one peak
        [~,ipk] = max(PKS);
        pk_height(iUn,ist) = PKS(ipk);
        pk_time(iUn,ist) = LOCS(ipk);
        
        % Find time before surpassing halfMax, and time after falling below        
        for ims = [LOCS(ipk):-1:1 pdms:-1:LOCS(ipk)]
            if ThisData(iUn,ims,ist)<halfMax
                t_onset(iUn,ist) = ims;
                break
            end
        end
        
        for ims = [LOCS(ipk):pdms 1:LOCS(ipk)]
            if ThisData(iUn,ims,ist)<halfMax
                t_offset(iUn,ist) = ims;
                break
            end
        end
        
        % Calculate peak duration
        if t_offset(iUn,ist)>=t_onset(iUn,ist)
            t_up(iUn,ist) = t_offset(iUn,ist)-t_onset(iUn,ist);
        else
            t_up(iUn,ist) = t_offset(iUn,ist)+pdms-t_onset(iUn,ist);
        end
        
        % Time of midpoint of peak
        t_mid(iUn,ist) = t_onset(iUn,ist)+t_up(iUn,ist)/2;
        if t_mid(iUn,ist)>pdms
            t_mid(iUn,ist) = t_mid(iUn,ist)-pdms;
        end
        
        
        % Add to plot
        subplot(2,3,ist); hold on
%         plot(pk_height(iUn,ist), t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',10)
%         plot(t_onset(iUn,ist), t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',pk_height(iUn,ist))
        plot(t_mid(iUn,ist),t_up(iUn,ist),'.','Color',colors(ist,:),'MarkerSize',10*(1+log10(pk_height(iUn,ist))) )
        
%         plot(t_onset(iUn,ist)+[0 t_up(iUn,ist)],[halfMax halfMax],'Color',colors(ist,:),'LineWidth',6)
        
    end %ist
end %iUn

% Finish plots
for ist = 1:size(ThisData,3)
    pdms = SegDurs(ist);
    
    subplot(2,3,ist); hold on
    axis square
    xlim([0 pdms])
    ylim([0 pdms])
    set(gca,'xtick',[0 pdms],'ytick',[0 pdms])
    xlabel('Mid time (ms)')
    ylabel('Peak duration (ms)')
    title(kfns{ist})
end

print_eps_kp(gcf,fullfile(savedir,'LatDur_Mids'))


ic = sum(isnan(t_up),2)>1;
t_up(ic,:) = nan;
avgPropUp = mean(t_up./repmat(SegDurs',[size(t_up,1) 1]),2);
figure; 
histogram(avgPropUp,0:0.02:1)

save(fullfile(fn.figs,'ClassContext','Speech','RawData','avgPropUp'),'avgPropUp','-v7.3')

