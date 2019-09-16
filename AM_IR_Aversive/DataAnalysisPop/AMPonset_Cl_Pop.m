function Data = AMPonset_Cl_Pop(JustPlot)
%
%  AMPonset_Cl_Pop
%   
%   By getting the raster data from MatchedMPs, cannot extract within
%   session simultaneous trial info.
% 
%  KP, 2019-09


global fn AMrates trMin Iterations

%!!!!!!!!!!!!!!!!!!!
Iterations  =  50;
%!!!!!!!!!!!!!!!!!!!
% iSess       =  5;
%!!!!!!!!!!!!!!!!!!!

if nargin<1
    JustPlot = 0;
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


%% Select Pdc MPs to use
% Find the period that occurs latest in the Pdc stimulus. Save the MPHrate
% table indices to make sure that it's the same period for every unit.

PdcMPindices = nan(5,2);

for irate = 1:5
    PdcMPindices(irate,:) = Data(21,irate).data(2,1).indices;
%     Data(21,irate).data(4,1)
end


%% Select which units to plot

% Index units by # SUs per session
% [Sessions,~,idxSess] = unique(UnitInfo(:,1:2),'stable');
% y=histogram(idxSess);
% [Nuns,iss] = sort(y.Values,'descend');
% SortedSessions = Sessions(iss,:);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Which units to collect data from?
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % all
theseUnits = 0:size(Data,1); 
% 
% % single session
% theseUnits = find(strcmp(UnitInfo.Subject,SortedSessions.Subject{iSess}) & strcmp(UnitInfo.Session,SortedSessions.Session{iSess}))';
% append_str = num2str(theseUnits(1));

% just RS
theseUnits = find([UnitInfo.TroughPeak]>=0.5)';
append_str = 'RS';
% 
% % just NS
% theseUnits = find([UnitInfo.TroughPeak]<0.5)';
% append_str = 'NS';


%%
if JustPlot
    % Load population d' results
%     load(fullfile(fn.processed,'MPcontext','Pop','Pop_dPrimes'))
else
    keyboard
    % write classifier here
    runPopClass_context
end


% Figure settings
set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');   %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Savedir for figure
savedir = fullfile(fn.processed,'AMPonset','Pop');
if ~exist(savedir,'dir')
    mkdir(savedir)
end



%% Compare population decoding to individual units

all_dPrimes = [];

hf=figure;
hold on
set(hf,'Position',fullscreen)

for irate = 1:size(Data,2)
    each_dPrime = [];
    
    for iUn = theseUnits
        
        if iUn==0 || isempty(Data(iUn,irate).data) || ~isfield(Data(iUn,irate).data,'IRid')
            continue
        end
        for ip = 1:size(Data(iUn,irate).data,1)
            each_dPrime = [each_dPrime; Res(iUn,irate).L1o.dprime(ip,2)];
            all_dPrimes = [all_dPrimes; Res(iUn,irate).L1o.dprime(ip,2)];
        end %ip
    end %iUn
    
    subplot(2,3,irate); hold on
    ndp = sum(~isnan(each_dPrime));
    plot([0 ndp+1],[0 0],':k')
    plot((1:numel(Pop_dPrimes(:,irate)))+ndp/2,Pop_dPrimes(:,irate),'om','MarkerSize',10)
    plot([0 ndp+1],[mean(Pop_dPrimes(:,irate)) mean(Pop_dPrimes(:,irate))],'-m','LineWidth',2)
    plot(sort(each_dPrime),'.k')
    axis square
    xlim([0 ndp+1])
    ylim([-3 6])
    title([num2str(AMrates(irate)) ' Hz'])
    if irate==1
        xlabel('MP comarison #')
        ylabel('dprime')
    end
end %irate


% All rates together
subplot(2,3,6); hold on
ndp = sum(~isnan(all_dPrimes));
plot([0 ndp+1],[0 0],':k')
hold on
plot((1:numel(Pop_dPrimes))+ndp/2,Pop_dPrimes(:),'ob','MarkerSize',10)
plot([0 ndp+1],[mean(Pop_dPrimes(:)) mean(Pop_dPrimes(:))],'-b','LineWidth',2)
plot(sort(all_dPrimes),'.k')
axis square
xlim([0 ndp+1])
ylim([-3 6])
title('All rates')

suptitle('MP context discrimination, SU vs Pop')

% Save
print_eps_kp(hf,fullfile(savedir,['SU-vs-Pop_' append_str]))



end


