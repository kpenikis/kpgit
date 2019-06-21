function pc_FRdependence
% 

global AMrates rateVec_AC rateVec_DB trMin 
close all

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
halfscreen = [1 scrsz(4) scrsz(3)/2 scrsz(4)];


% Colors
colors = [  84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;

% Other settings
histbinsize = 0.025;
trMin       =  10;

AMrates     = [2 4 8 16 32];

q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;


%%  FR HISTORY DEPENDENCE FOR HIGH D' 

dp_thresh   = 1;
alpha_thr   = 0.01;
thisMet     = 'Prev100msFR';
nIterations = 1000;

hf=figure; 
set(hf,'Position',halfscreen,'NextPlot','add')
% plot([10^-2 10^2],[10^-2 10^2],'-c')
% set(gca,'xscale','log','yscale','log','xlim',[10^-2 10^2],'ylim',[10^-2 10^2])
plot([-5 50],[-5 50],'-r')
set(gca,'xscale','linear','yscale','linear','xlim',[-5 50],'ylim',[-5 50])
axis square


%% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
%            Shuffled spikes
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

pFRsh_sig  = 0;

for iUn = 1:size(Data,1)
    %-----------------------------------
    % Get trial-by-trial firing history -- IF hasn't been saved that way already
    if length(Data(iUn,2).data(1).(thisMet))<2
        keyboard
        DataSh(iUn,:) = getFRhist( DataSh(iUn,:), UnitData(iUn), spkshift );
    end
    %-----------------------------------
    for irate = 1:size(Data,2)
    
        if isempty(Data(iUn,irate).data), continue, end
        
        dprimes = Data(iUn,irate).Res_L1o.dprime(:,2);
        high_dp = find(dprimes>dp_thresh)';
    
        for idp = high_dp
            
            % Make sure Data struct updated
            if length(Data(iUn,irate).data(1).(thisMet))<2
                keyboard
            end
            
            clear p r rp
                        
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);     
            
            % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
            %  Shuffled spikes
            % ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
            
            ntp = size(DataSh(iUn,irate).data(1).(thisMet),1);
            nto = size(DataSh(iUn,irate).data(id).(thisMet),1);
            
            % Bootstrap to compare distributions
            p=ones(nIterations,1);
            dm=ones(nIterations,1);
            
            if ntp>nto
                bs_data = DataSh(iUn,irate).data(1).(thisMet);
                ss_data = DataSh(iUn,irate).data(id).(thisMet);
                nt = nto;
            elseif ntp<nto
                bs_data = DataSh(iUn,irate).data(id).(thisMet);
                ss_data = DataSh(iUn,irate).data(1).(thisMet);
                nt = ntp;
            end
            
            for ii = 1:nIterations
                p(ii)  = ranksum(ss_data,bs_data(randi(length(bs_data),1,nt)));
                dm(ii) = median(ss_data) - median(bs_data(randi(length(bs_data),1,nt)));
            end
            qr = quantile(p,0.025);
            qd = quantile(dm,[0.025 0.975]);
            if qd(1)>0 || qd(2)<0    %qr(1)<=alpha_thr
                sigFR = 1;
            else
                sigFR = 0;
            end
            
            
            % Plot results of same datapoints but shuffled spikes
            x_mean   = mean(DataSh(iUn,irate).data(1).(thisMet));
            y_mean   = mean(DataSh(iUn,irate).data(id).(thisMet));
            x_std    = std(DataSh(iUn,irate).data(1).(thisMet));
            y_std    = std(DataSh(iUn,irate).data(id).(thisMet));
            
            figure(hf); hold on
            
            if sigFR
                plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-','LineWidth',1,'Color','k')
                plot([x_mean x_mean],y_mean + y_std*[-1 1],'-','LineWidth',1,'Color','k')
                plot(x_mean,y_mean,'.','MarkerSize',50,'Color','k')
                
                pFRsh_sig = pFRsh_sig+1;
            else
                plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-','LineWidth',1,'Color',0.7*[1 1 1])
                plot([x_mean x_mean],y_mean + y_std*[-1 1],'-','LineWidth',1,'Color',0.7*[1 1 1])
                plot(x_mean,y_mean,'.','MarkerSize',50,'Color',0.7*[1 1 1])
            end
                        
        end %high_dp
    end %irate
end %iUn


%% ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
%            Real spikes
%  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

Nh         = 0;
Nl         = 0;
pFR_sig    = 0;
rp_Dh_s      = [];
rp_Dh_ns      = [];

for iUn = 1:size(Data,1)
    %-----------------------------------
    % Get trial-by-trial firing history -- IF hasn't been saved that way already
    if length(Data(iUn,2).data(1).(thisMet))<2
        keyboard
        Data(iUn,:)   = getFRhist( Data(iUn,:),   UnitData(iUn), spkshift );
    end
    %-----------------------------------
    for irate = 1:size(Data,2)
    
        if isempty(Data(iUn,irate).data), continue, end
        
        dprimes = Data(iUn,irate).Res_L1o.dprime(:,2);
        high_dp = find(dprimes>dp_thresh)';
        low_dp  = find(dprimes<=dp_thresh)';
        
        for idp = low_dp
            
            Nl = Nl+1;
            
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
%             plot(DataSh(iUn,irate).data(1).(thisMet),DataSh(iUn,irate).data(id).(thisMet),'o','Color',0.5*[1 1 1],'LineWidth',2)
            
        end %low_dp
    
        for idp = high_dp
            
            % Make sure Data struct updated
            if length(Data(iUn,irate).data(1).(thisMet))<2
                keyboard
            end
            
            clear p r rp
            
            Nh = Nh+1;
            
            id = Data(iUn,irate).Res_L1o.dprime(idp,1);     
            
            
            % - - - - - - - - - -
            %    Real spikes
            % - - - - - - - - - -
            
            ntp = size(Data(iUn,irate).data(1).(thisMet),1);
            nto = size(Data(iUn,irate).data(id).(thisMet),1);
            
            % Bootstrap to compare distributions
            p=ones(nIterations,1);
            
            if ntp>nto
                bs_data = Data(iUn,irate).data(1).(thisMet);
                ss_data = Data(iUn,irate).data(id).(thisMet);
                nt = nto;
            elseif ntp<nto
                bs_data = Data(iUn,irate).data(id).(thisMet);
                ss_data = Data(iUn,irate).data(1).(thisMet);
                nt = ntp;
            end
            
            for ii = 1:nIterations
                p(ii)  = ranksum(ss_data,bs_data(randi(length(bs_data),1,nt)));
                dm(ii) = median(ss_data) - median(bs_data(randi(length(bs_data),1,nt)));
            end
            qr = quantile(p,0.025);
            qd = quantile(dm,[0.025 0.975]);
            if qd(1)>0 || qd(2)<0    %qr(1)<=alpha_thr
                sigFR = 1;
            else
                sigFR = 0;
            end
            
            
            % Plot results of datapoints with high dprime
            x_mean   = mean(Data(iUn,irate).data(1).(thisMet));
            y_mean   = mean(Data(iUn,irate).data(id).(thisMet));
            x_std    = std(Data(iUn,irate).data(1).(thisMet));
            y_std    = std(Data(iUn,irate).data(id).(thisMet));
            
            figure(hf); hold on
            
            if sigFR
                
                plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-','LineWidth',1,'Color','b')
                plot([x_mean x_mean],y_mean + y_std*[-1 1],'-','LineWidth',1,'Color','b')
                plot(x_mean,y_mean,'.','MarkerSize',50,'Color','b') 
                
                Data(iUn,irate).data(id).PrecFR_sig = 1;
                pFR_sig = pFR_sig+1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%   Plot trial by trial relationship for sig units  %%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                hfTr(pFR_sig) = figure; 
                plot([0 50],[0 50],'-k')
                hold on
                ip=plot(Data(iUn,irate).data(1).(thisMet), sum(Data(iUn,irate).data(1).raster,2)/(1000/Data(iUn,irate).AMrate)*1000,...
                    '.','MarkerSize',30,'Color',colors(irate,:));
                plot(Data(iUn,irate).data(id).(thisMet), sum(Data(iUn,irate).data(id).raster,2)/(1000/Data(iUn,irate).AMrate)*1000,...
                    'o','MarkerSize',15,'LineWidth',2,'Color',colors(irate,:))
                axis square
                title(sprintf('Trial FR effects for  iUn: %i  irate: %i  id: %i',iUn,irate,id))
                
                % Add correlation of Preceding FR to current response
                [r,rp] = corrcoef([Data(iUn,irate).data(1).(thisMet); Data(iUn,irate).data(id).(thisMet)],...
                    [sum(Data(iUn,irate).data(1).raster,2)/(1000/Data(iUn,irate).AMrate)*1000; sum(Data(iUn,irate).data(id).raster,2)/(1000/Data(iUn,irate).AMrate)*1000]);
                legend(ip,sprintf('r = %0.3f\np = %0.2g',r(1,2),rp(1,2)))
                
                % Store p value of correlation
                rp_Dh_s = [rp_Dh_s rp(1,2)];
                
                
            else % non-significant preceding rate distributions
                
                plot(x_mean + x_std*[-1 1],[y_mean y_mean],'-','LineWidth',1,'Color',[0.65 0.82 1])
                plot([x_mean x_mean],y_mean + y_std*[-1 1],'-','LineWidth',1,'Color',[0.65 0.82 1])
                plot(x_mean,y_mean,'.','MarkerSize',50,'Color',[0.65 0.82 1])
                
                % Calculate correlation of Preceding FR to current response
                [~,rp] = corrcoef([Data(iUn,irate).data(1).(thisMet); Data(iUn,irate).data(id).(thisMet)],...
                    [sum(Data(iUn,irate).data(1).raster,2)/(1000/Data(iUn,irate).AMrate)*1000; sum(Data(iUn,irate).data(id).raster,2)/(1000/Data(iUn,irate).AMrate)*1000]);                
                % Store p value of correlation
                rp_Dh_ns = [rp_Dh_ns rp(1,2)];
                
            end
                   
        end %high_dp
    end %irate
end %iUn


%% Finish and save plot of avg preceding rates (high dprime cases only)

figure(hf); hold on
xlabel(sprintf('%s for Pdc MPH\ngrey: shuffled spikes',thisMet))
ylabel(sprintf('%s for other MPH',thisMet))
title(['Change in response as a function of ' thisMet ', N=' num2str(Nh) ', dp thresh=' num2str(dp_thresh)])

print_eps_kp(hf,fullfile(fn.figs,'MPHclass',['FRdependence_dp' num2str(dp_thresh*10) '_' thisMet] ))


%% Print all results to a txt file

filename = fullfile(fn.figs,'MPHclass',sprintf('MPHclassResults_dp%i_%s.txt',dp_thresh*10,thisMet));

fileID = fopen(filename,'w');

fprintf(fileID,'##MPH CLASSIFIER RESULTS \n\n  %i/%i MPH comparisons have d-primes > %0.1f \n', Nh,Nh+Nl,dp_thresh );
fprintf(fileID,'  Of high d-primes, %i/%i show sig diff distributions of %s \n',pFR_sig,Nh,thisMet );
fprintf(fileID,'  (shuffled spikes) %i/%i \n',pFRsh_sig,Nh );

% Compare high to low dprime cases for trial by trial FR history relationship 
fprintf(fileID,'  Of cases with sig %s, %i/%i have a sig correlation between %s and MPH resp on a trial by trial basis \n',...
    thisMet, sum(rp_Dh_s<alpha_thr), pFR_sig, thisMet );
fprintf(fileID,'  Of cases with NON sig %s, %i/%i have a sig correlation between %s and MPH resp on a trial by trial basis \n\n',...
    thisMet, sum(rp_Dh_ns<alpha_thr), Nh-pFR_sig, thisMet );

fclose(fileID);

type(filename)


end


