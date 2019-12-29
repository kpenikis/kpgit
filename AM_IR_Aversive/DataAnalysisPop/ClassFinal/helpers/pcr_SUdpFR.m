
%% ~~~~ Distribution of PC vals for each stimulus, by mean FR

nSpks  = (mean(sum(mean(CTTS,3,'omitnan'),2,'omitnan'),4));
nSpkRS = (mean(sum(mean(CTTS(iRS,:,:,:),3,'omitnan'),2,'omitnan'),4));
nSpkNS = (mean(sum(mean(CTTS(iNS,:,:,:),3,'omitnan'),2,'omitnan'),4));

figure;
subplot(1,4,1)
hold on


% Manually make boxplots
q5  = quantile(CReach.dprime,0.05);
q25 = quantile(CReach.dprime,0.25);
q75 = quantile(CReach.dprime,0.75);
q95 = quantile(CReach.dprime,0.95);

plot([1 1],[q5 q95],'-','Color',[0.5 0.7 1],'LineWidth',6)
fill(1+[0.3 0.3 -0.3 -0.3],[q75 q25 q25 q75],[0.5 0.7 1],'EdgeColor','none')

ylabel('dprime')
set(gca,'Color','none')
box off
ylim([-0.5 2.5])
plotSpread(CReach.dprime,'distributionIdx',ones(size(CReach.dprime)),'distributionColors','k','showMM',3)


subplot(1,4,2:4)
%         plot(nSpks,CReach.dprime,'k.')
hold on
plot(nSpkRS,CReach.dprime(iRS),'r.')
plot(nSpkNS,CReach.dprime(iNS),'b.')
ylim([-0.5 2.5])
set(gca,'ytick',[],'Color','none')
box off
xlabel('N spks')
%         title('SU dprimes vs mean N spks per stim')


% Stats
[r_RS,p_RS] = corr(nSpkRS,CReach.dprime(iRS),'Type','Spearman');
[r_NS,p_NS] = corr(nSpkNS,CReach.dprime(iNS),'Type','Spearman');
[r,p]       = corr(nSpks,CReach.dprime,'Type','Spearman');

text(40,-0.3,sprintf('Spearman r=%0.2f, p=%0.3e',r,p))


savename = sprintf('SU_dps_%s_%s',varPar,whichStim);

keyboard

print_eps_kp(gcf,fullfile(savedir,savename))

