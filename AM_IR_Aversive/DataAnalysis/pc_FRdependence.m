function pc_FRdependence
% 

global AMrates rateVec_AC rateVec_DB trMin 

%% 
% Load Unit files
fn = set_paths_directories;
q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------
[sigUnits,UnitData] = identifyResponsiveUnits(UnitData);

% Load classifier results
load(fullfile(fn.processed,'MPHclassifier','ClassData'));
q=load(fullfile(fn.processed,'MPHclassifier','ClassData_shuff'));
DataSh = q.Data;
clear q

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)
scrsz = get(0,'ScreenSize');  %[left bottom width height]
tallrect   = [1 scrsz(4) scrsz(3)/4 scrsz(4)];
longrect   = [1 scrsz(4) scrsz(3) scrsz(4)/4];
halfrect   = [1 scrsz(4) scrsz(3) scrsz(4)/2];
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];


% Other settings
histbinsize = 0.025;
trMin       =  10;

AMrates     = [2 4 8 16 32];

q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;


%%  FR HISTORY DEPENDENCE FOR HIGH D' 

dp_thresh = 1.5;
thisMet = 'Prev100msFR';

hf=figure; 
set(hf,'Position',fullscreen,'NextPlot','add')
% plot([10^-2 10^2],[10^-2 10^2],'-c')
% set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
plot([-5 100],[-5 100],'-r')
set(gca,'xscale','linear','yscale','linear','xlim',[-5 100],'ylim',[-5 100])
axis square

Nh = 0;
Nl = 0;
for iUn = 1:size(Data,1)
    %-----------------------------------
    % Get trial-by-trial firing history -- IF hasn't been saved that way already
    if length(Data(iUn,2).data(1).(thisMet))<2
        keyboard
        Data(iUn,:)   = getFRhist( Data(iUn,:),   UnitData(iUn), spkshift );
        DataSh(iUn,:) = getFRhist( DataSh(iUn,:), UnitData(iUn), spkshift );
    end
    %-----------------------------------
    for irate = 1:size(Data,2)
    
        if isempty(Data(iUn,irate).data), continue, end
        
        dprimes = Data(iUn,irate).Res_L1o.dprime(:,2);
        high_dp = find(dprimes>dp_thresh)';
        low_dp  = find(dprimes<=dp_thresh)';
        
%         for idp = low_dp
%             
%             Nl = Nl+1;
%             
%             id = Data(iUn,irate).Res_L1o.dprime(idp,1);            
%             
%             x_mean   = mean(Data(iUn,irate).data(1).(thisMet));
%             y_mean   = mean(Data(iUn,irate).data(id).(thisMet));
%             x_std    = std(Data(iUn,irate).data(1).(thisMet));
%             y_std    = std(Data(iUn,irate).data(id).(thisMet));
%             
%             figure(hf); hold on
%             plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-','Color',0.7*[1 1 1],'LineWidth',2)
%             plot([x_mean x_mean],y_mean + y_std*[-1 1],'-','Color',0.7*[1 1 1],'LineWidth',2)
%             plot(x_mean,y_mean,'o','Color',0.7*[1 1 1],'LineWidth',2)
% %             plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'o','Color',0.5*[1 1 1],'LineWidth',2)
%             
%         end %low_dp
        for idp = high_dp
            % Make sure Data struct updated
            if length(Data(iUn,irate).data(1).(thisMet))<2
                keyboard
            end
            Nh = Nh+1;
            
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);            
            
            % Plot results of datapoints with high dprime
            x_mean   = mean(Data(iUn,irate).data(1).(thisMet));
            y_mean   = mean(Data(iUn,irate).data(id).(thisMet));
            x_std    = std(Data(iUn,irate).data(1).(thisMet));
            y_std    = std(Data(iUn,irate).data(id).(thisMet));
            
            figure(hf); hold on
            plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-b','LineWidth',2)
            plot([x_mean x_mean],y_mean + y_std*[-1 1],'-b','LineWidth',2)
            plot(x_mean,y_mean,'ob','LineWidth',2)
            
            % Plot results of same datapoints but shuffled spikes
            x_mean   = mean(DataSh(iUn,irate).data(1).(thisMet));
            y_mean   = mean(DataSh(iUn,irate).data(id).(thisMet));
            x_std    = std(DataSh(iUn,irate).data(1).(thisMet));
            y_std    = std(DataSh(iUn,irate).data(id).(thisMet));
            
            plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-k','LineWidth',2)
            plot([x_mean x_mean],y_mean + y_std*[-1 1],'-k','LineWidth',2)
            plot(x_mean,y_mean,'ok','LineWidth',2)
            
%             plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'ok','LineWidth',2)
            
        end %high_dp
        
    end %irate
end %iUn


% Finish and save fig 1
figure(hf); hold on
title(['Change in response as a function of ' thisMet ', N=' num2str(Nh) ', dp thresh=' num2str(dp_thresh)])
set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
axis square
plot([10^-2 10^2],[10^-2 10^2],'-c')

% print_eps_kp(hf,fullfile(fn.figs,'MPHclass',['RespDependence_dpAbove_' num2str(dp_thresh*10) '_' thisMet '_PreToPre']))



% Resave Data structure !!



% Stat test for what datapoints are off unity line
% Take these dps, look at trial by trial FR history relationship 





end


