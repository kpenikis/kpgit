function Data = PCl_MPcontext(JustPlot)
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

% SU classifier results
q = load(fullfile(fn.processed,'MPHclassifier','MPcontextSUclassRes'));
Res = q.Res;
clear q


%% Select which units to plot

% Index units by # SUs per session
[Sessions,~,idxSess] = unique(UnitInfo(:,1:2),'stable');
y=histogram(idxSess);
[Nuns,iss] = sort(y.Values,'descend');
SortedSessions = Sessions(iss,:);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Which units to collect data from?
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%   ALL
% theseUnits = 0:size(Data,1); 
% append_str = '';

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
k          = setSUdpCutoff(Res,0.95);
theseUnits = 1:size(Data,1);
append_str = 'SUabove_q95';


%%
if JustPlot
    % Load population d' results
    load(fullfile(fn.processed,'MPcontext','Pop','Pop_dPrimes'))
else
    fprintf('Running classifier on population data\n')
    runPopClass_context
end


% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Savedir for figure
savedir = fullfile(fn.processed,'MPcontext','Pop');
if ~exist(savedir,'dir')
    mkdir(savedir)
end


%% Compare population decoding to individual units


% Get full poulation decoding performance for comparison
if strncmp(append_str,'exclSU',6) || strncmp(append_str,'SUabov',6)
    p=load(fullfile(fn.processed,'MPcontext','Pop','Pop_dPrimes'));
end


all_dPrimes    = [];
all_dPrimes_ex = [];

hf=figure;
hold on
set(hf,'Position',fullscreen)

for irate = 1:size(Data,2)
    each_dPrime = [];
    excl_dPrime = [];
    
    for iUn = theseUnits
        
        if iUn==0 || isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'IRid')
            continue
        end
        for ip = 1:size(Data(iUn,irate).data,1)
            
            if exist('q','var') && Res(iUn,irate).L1o.dprime(ip,2)>q(irate)
                excl_dPrime = [excl_dPrime; Res(iUn,irate).L1o.dprime(ip,2)];
            elseif exist('k','var') && Res(iUn,irate).L1o.dprime(ip,2)<k(irate)
                excl_dPrime = [excl_dPrime; Res(iUn,irate).L1o.dprime(ip,2)];
            elseif ~exist('k','var') || ( exist('k','var') && Res(iUn,irate).L1o.dprime(ip,2)>k(irate) )
                each_dPrime = [each_dPrime; Res(iUn,irate).L1o.dprime(ip,2)];
            end
            
        end %ip
    end %iUn
    
    all_dPrimes    = [all_dPrimes; each_dPrime];
    all_dPrimes_ex = [all_dPrimes_ex; excl_dPrime];
    
    subplot(2,3,irate); hold on
    ndp = sum(~isnan(each_dPrime));
%     plot([0 ndp+1+sum(~isnan(all_dPrimes_ex))],[0 0],'k','LineWidth',0.5)
    plot((1:numel(Pop_dPrimes(:,irate)))+ndp/2,Pop_dPrimes(:,irate),'om','MarkerSize',10)
    plot([0 ndp+1],[mean(Pop_dPrimes(:,irate)) mean(Pop_dPrimes(:,irate))],'-m','LineWidth',2)
    plot(sort(each_dPrime),'.k')
    plot(ndp+(1:length(excl_dPrime)),sort(excl_dPrime),'.y')
    if exist('p','var')
        plot([0 ndp+1+sum(~isnan(excl_dPrime))],[mean(p.Pop_dPrimes(:,irate)) mean(p.Pop_dPrimes(:,irate))],':m','LineWidth',1)
    end
    axis square
    xlim([0 ndp+1+sum(~isnan(excl_dPrime))])
    ylim([-1 6])
    title([num2str(AMrates(irate)) ' Hz'])
    if irate==1
        xlabel('MP comparison #')
        ylabel('dprime')
    end
end %irate


% All rates together
subplot(2,3,6); hold on
ndp = sum(~isnan(all_dPrimes));
% plot([0 ndp+1+sum(~isnan(all_dPrimes_ex))],[0 0],'k','LineWidth',0.5)
hold on
plot((1:numel(Pop_dPrimes))+ndp/2,Pop_dPrimes(:),'ob','MarkerSize',10)
plot([0 ndp+1],[mean(Pop_dPrimes(:)) mean(Pop_dPrimes(:))],'-b','LineWidth',2)
plot(sort(all_dPrimes),'.k')
plot(ndp+(1:length(all_dPrimes_ex)),sort(all_dPrimes_ex),'.y')
if exist('p','var')
    plot([0 ndp+1+sum(~isnan(all_dPrimes_ex))],[mean(p.Pop_dPrimes(:)) mean(p.Pop_dPrimes(:))],':b','LineWidth',1)
end
axis square
xlim([0 ndp+1+sum(~isnan(all_dPrimes_ex))])
ylim([-1 6])
title('All rates')

suptitle('MP context discrimination, SU vs Pop')

beep

%% Save
print_eps_kp(hf,fullfile(savedir,['SU-vs-Pop_' append_str]))



end


