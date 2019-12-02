function grpClassResults


TempSize  = 19;
binsize   = 10;

fn = set_paths_directories('','',1);

savedir = fullfile(fn.figs,'StimClass','Grouped');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

q=load(fullfile(savedir,sprintf('PCMat_grp_random_%ims_tempsz%i',binsize,TempSize)));
PCrand   = q.RandGrpPCMat;

GroupStr = {'phase' 'peakFR' 'tuning' 'dynRange'};
PCgrp    = nan(size(PCrand,1),size(PCrand,1),numel(GroupStr));
for ig = 1:numel(GroupStr)
    q=load(fullfile(savedir,sprintf('PCMat_grp_%s_%ims_tempsz%i',GroupStr{ig},binsize,TempSize)));
    PCgrp(:,:,ig) = q.PCMat;
end


% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');     %[left bottom width height]
widescreen  = [1 scrsz(4)/2 scrsz(3) scrsz(4)/2];

colgrps = cmocean('phase',5);
grpjit  = -0.15:0.1:0.15;


%% Make figure

hf=figure; 
set(hf,'Position',widescreen)
subplot(1,5,1:4); 
hold on

% PERCENT CORRECT

for ist = 1:size(PCrand,1)
    
%     boxplot(permute(PCrand(ist,ist,:),[3 1 2]),ist*ones(size(PCrand,3),1))
    plotSpread(permute(PCrand(ist,ist,:),[3 1 2]),'distributionIdx',ist*ones(size(PCrand,3),1),'distributionColors','k')
    
    clear ip
    for ig = 1:numel(GroupStr)
        ip(ig)=plot(ist+grpjit(ig),PCgrp(ist,ist,ig),'x','Color',colgrps(ig,:),...
            'LineWidth',3,'MarkerSize',20);
    end
    
end
legend(ip,GroupStr,'Location','southeast','Color','none')
set(gca,'xticklabel',[{'Warn'} {'2Hz'} {'4Hz'} {'8Hz'} {'16Hz'} {'32Hz'}],...
    'xtick',1:ist,'tickdir','out','Color','none')

ylabel('Proportion Correct')
ylim([0 1])
% title(sprintf('Proportion Correct by stimulus: euclidean classifier, leave 1 out, conv %ims',binsize))
title('Proportion Correct by stimulus')


% DPRIME

dprimes_rand = nan(size(PCrand,3),1);
for ish = 1:size(PCrand,3)
    thisMat = PCrand(:,:,ish);
    dprimes_rand(ish) = norminv(mean(diag(thisMat)),0,1) - norminv(mean(thisMat(~diag(diag(thisMat)))),0,1); 
end

subplot(1,5,5)
hold on
plotSpread(dprimes_rand,'distributionColors','k')

dprime_grps = nan(size(PCgrp,3),1);
for ig = 1:numel(GroupStr)
    dprime_grps(ig) = norminv(mean(diag(PCgrp(:,:,ig))),0,1) - norminv(mean(PCgrp(~diag(diag(PCgrp(:,:,ig))))),0,1); 
    ip(ig)=plot(1+ig,dprime_grps(ig),'x','Color',colgrps(ig,:),...
        'LineWidth',3,'MarkerSize',20);
end
% legend(ip,GroupStr,'Location','southeast','Color','none')
set(gca,'xticklabel',[],'xtick',1:(ig+1),'tickdir','out','Color','none')
ylabel('d-prime')
ylim([0 3])
xlim([0 ig+2])
title('d-primes')
xlabel('Grouping type')



print_eps_kp(hf,fullfile(savedir,sprintf('Results_temp%i_bin%i',TempSize,binsize)))



end %function


