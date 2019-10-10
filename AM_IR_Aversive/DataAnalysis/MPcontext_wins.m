function Data = MPcontext_wins(JustPlot)
%
%  PCl_MPcontext(JustPlot)
%   Compare discriminability of MP context by Population and SUs. 
%   
%   opt argin: JustPlot
%       if JustPlot=0, re-runs classifier on population data.
%   
%   Must select "theseUnits" to compare.
%
%  KP, 2019-09
%

global fn AMrates trMin Iterations

%!!!!!!!!!!!!!!!!!!!
Iterations  =  5000;
%!!!!!!!!!!!!!!!!!!!

if nargin<1
    JustPlot = 1;
end

%% Set up

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load matched MP data
q = load(fullfile(fn.processed,'MatchedMPs','MatchedMPdata'));
Data = q.Data;
clear q

% SU classifier results  --  whole period
q = load(fullfile(fn.processed,'MPHclassifier','MPcontextSUclassRes'));
Res_pd = q.Res;
clear q

% SU classifier results  --  62 ms
q = load(fullfile(fn.processed,'MPcontext','MPcontextSU_62'));
Res_62 = q.Res;
clear q


%% Select which units to plot

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Which units to collect data from?
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   ALL
theseUnits = 0:size(Res_62,1); 
append_str = '';

%   Single session
% iSess       =  5;
% theseUnits = find(strcmp(UnitInfo.Subject,SortedSessions.Subject{iSess}) & strcmp(UnitInfo.Session,SortedSessions.Session{iSess}))';
% append_str = num2str(theseUnits(1));

%   just RS
% theseUnits = find([UnitInfo.TroughPeak]>=0.5)';
% append_str = 'RS';

%   just NS
% theseUnits = find([UnitInfo.TroughPeak]<0.5)';
% append_str = 'NS';

%   Remove best SU 
% q          = setSUdpCutoff(Res,0.85);
% theseUnits = 1:size(Data,1);
% append_str = 'exclSU_q85';

%   Only keep best SU 
% k          = setSUdpCutoff(Res,0.95);
% theseUnits = 1:size(Data,1);
% append_str = 'SUabove_q95';


%%

% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Savedir for figure
savedir = fullfile(fn.processed,'MPcontext');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Compare population decoding to individual units

all_dp_pd  = [];
all_dp_62  = [];

hf1=figure;
hold on
set(hf1,'Position',fullscreen)

hf2=figure;
hold on
set(hf2,'Position',fullscreen)

for irate = 1:size(Res_62,2)
    dPrime_pd = [];
    dPrime_62 = [];
    
    for iUn = theseUnits
        
        if iUn==0 || isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'IRid')
            continue
        end
        for ip = 1:size(Data(iUn,irate).data,1)
            
            dPrime_pd = [dPrime_pd; Res_pd(iUn,irate).L1o.dprime(ip,2)];
            dPrime_62 = [dPrime_62; Res_62(iUn,irate).L1o.dprime(ip,2)];
            
        end %ip
    end %iUn
    
    all_dp_pd    = [all_dp_pd; dPrime_pd];
    all_dp_62    = [all_dp_62; dPrime_62];
    
    [~,isorted] = sort(dPrime_pd);
    
    figure(hf1); hold on
    subplot(2,2,irate); hold on
    plot([-1 6],[-1 6],'k')
    plot(dPrime_pd(isorted),'.k')
    plot(dPrime_62(isorted),'.b')
    axis square
    xlim([0 length(dPrime_pd)+1])
    ylim([-1 6])
    title([num2str(AMrates(irate)) ' Hz'])
    if irate==3
        ylabel(sprintf('dprimes\nblack: whole period, blue: first 62 ms'))
        xlabel('MP comparison ranked')
    end
    
    figure(hf2); hold on
    plotSpread([dPrime_pd; dPrime_62],'distributionIdx',[ones(size(dPrime_pd)); 2.*ones(size(dPrime_62))],'distributionColors',{'k' 'b'},'xValues',irate+[-0.2 0.2])
    
end %irate

set(gca,'xtick',1:4,'xticklabel',AMrates(1:4))
ylabel('dprime')


%% Save

print_eps_kp(hf1,fullfile(savedir,'SU_fullPd-vs-62ms_ranked'))
print_eps_kp(hf2,fullfile(savedir,'SU_fullPd-vs-62ms_spread'))



end


