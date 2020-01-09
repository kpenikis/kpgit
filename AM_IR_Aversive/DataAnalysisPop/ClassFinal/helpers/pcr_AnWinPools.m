
close all

varPar       = 'AnDurPools';
whichCells   = 'Mid10RS'; % BestRS SkipBestRS 
whichStim    = 'Speech';

fn = set_paths_directories('','',1);
switch whichStim
    case {'AC' 'DB'}
        rootdir = fullfile(fn.figs,'ClassAM');
    case 'Speech'
        rootdir = fullfile(fn.figs,'ClassSpeech');
end
figsavedir = fullfile(rootdir,whichStim,varPar,whichCells);


% Get pooled cells results
if ~exist('CR','var')
    fprintf('no Results table in workspace, so loading it..\n')
    tablesavename = sprintf('CR_v%s_%s.mat',varPar,whichCells);
    q=load(fullfile(figsavedir,tablesavename));
    CR = q.CR;
    clear q
end

% Load Full results
tablesavename = sprintf('CR_v%s_%s.mat','Full','Mid10RS');
q=load(fullfile(rootdir,whichStim,'Full','Mid10RS',tablesavename));
CRfull = q.CR;
clear q

Pools = unique(CR.iC);
WinEnds = unique(CR.WinEnd);
%add whole 500 ms too (load another CR table) 

hf=figure; 
hold on

for ip = 1:numel(Pools)
    
    idx = find(CR.iC==Pools(ip));
    
    [WinEnds,isort] = sort(CR.WinEnd(idx)-500);
    
    % Add dp for full 500 ms duration
    dp_full = CRfull.dprime(CRfull.iC==Pools(ip));
    
    hp(ip)=plot ([WinEnds; 500] , [CR.dprime(idx(isort)); dp_full] ,'LineWidth',2);
    
end

legend(cellstr(num2str(Pools))','Location','best')
ylabel('dprime')
xlabel('An Win duration')
set(gca,'Color','none')

keyboard
print_eps_kp(gcf,fullfile(figsavedir,'dp_AnWinPools'))


