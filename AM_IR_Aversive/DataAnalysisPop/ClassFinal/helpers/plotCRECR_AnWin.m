% Plot neural discrimination vs envelope discrimination
%  as a function of increasing analysis window. 

q=load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/StimClassRcorr/CR_vAnDur_AC_all.mat');
CR = q.CR;
clear  q
q=load('/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/Env/EnvCR_vAnDur_AC_all.mat');
ECR = q.ECR;
clear q

theseStim = 1:8;
[WinEnds,iCR] = unique(CR.WinEnd);

StimPCs = nan(8,length(WinEnds));
StimDPs = nan(8,length(WinEnds));
for idur = 1:size(CR,1)
    for ist = theseStim
        
        data = mean(CR(idur,:).Results{:},3);
        StimPCs(ist,iCR(idur)) = data(ist,ist);
        StimDPs(ist,iCR(idur)) = norminv(data(ist,ist),0,1) - norminv(mean(data(ist,theseStim~=ist)),0,1);
        
    end
end

% ECR = ECR(cellfun(@length,ECR.Cells)==50,:);
[EWinEnds,iECR] = unique(ECR.WinEnd);
EnvPCs = nan(8,length(EWinEnds));
EnvDPs = nan(8,length(EWinEnds));
for idur = 1:size(ECR,1)
    for ist = theseStim
        
        data = mean(ECR(idur,:).Results{:},3);
        EnvPCs(ist,iECR(idur)) = data(ist,ist);
        EnvDPs(ist,iECR(idur)) = norminv(data(ist,ist),0,1) - norminv(mean(data(ist,theseStim~=ist)),0,1);
        
    end
end

% Set colors
colors = [ 150 150 150;...
            84  24  69;...
           120  10  41;...
           181   0  52;...
           255  87  51;...
           255 153   0]./255;
colors = [ colors; ...
            [37  84 156]./255 ;...
            [19 125 124]./255 ];

figure;
hold on
for ist = theseStim
    subplot(2,4,ist)
    hold on
    plot(WinEnds-500,StimPCs(ist,:),'LineWidth',4,'Color',colors(ist,:))
    plot(EWinEnds-500,EnvPCs(ist,:),':','LineWidth',2,'Color',colors(ist,:))
    ylim([0 1])
    set(gca,'Color','none')
    title(ist)
end

print_eps_kp(gcf,'/Volumes/GoogleDrive/My Drive/Sanes/DATADIR/AMaversive/Figures/ClassAM/CRECR_AnWin')



