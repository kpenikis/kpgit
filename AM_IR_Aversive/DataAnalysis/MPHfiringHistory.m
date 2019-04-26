function Data = classifyMPHcontext(USE_MEASURE)
%
%  classifyMPHcontext
%
%
%
%  KP, 2019-04
%


global fn AMrates rateVec_AC rateVec_DB trMin Iterations
fn = set_paths_directories([],[],1);

if nargin<1
    USE_MEASURE = 'FR'; 'spikes';
end

%!!!!!!!!!!!!!!!!!!!
trMin       =  10;
%!!!!!!!!!!!!!!!!!!!
Iterations  =  500;
%!!!!!!!!!!!!!!!!!!!
tVarBin     = 31;
%!!!!!!!!!!!!!!!!!!!
N=0;


%% Load Unit data files

fn = set_paths_directories('','',1);

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


%% Select a datapoint (or a few?) to highlight

ex_subj = 'xWWWf_253400';
ex_sess = 'GA';
ex_ch   = 16;
ex_clu  = 1;


%% Preallocate

pvalues = [];

%% Figure settings

set(0,'DefaultTextInterpreter','none')
set(0,'DefaultAxesFontSize',14)

scrsz = get(0,'ScreenSize');  %[left bottom width height]
fullscreen = [1 scrsz(4) scrsz(3) scrsz(4)];

% Set colors
colors = [ 0 200 150;...
    84  24  69;...
    120  10  41;...
    181   0  52;...
    255  87  51;...
    255 153   0]./255;
colors = [ colors; ...
    [37  84 156]./255 ;...
    [19 125 124]./255 ];


yvals = 2.^[0:7];


%% Step through Units

for iUn = 1:17%numel(UnitData)
        
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
    if UnitData(iUn).kw_p>0.05
        continue
    end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    subjcol     = [1 1 1];
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
    
    spkshift    = 0;%UnitData(iUn).IntTime_spk;
    
    
    % Load data files
    
%     if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info RateStream SoundStream SpoutStream
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
%     end
%     if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
%     end
    if ~isfield(Info,'artifact')
        continue
    end
    
    if ~exist('RateStream','var')
        keyboard
    end
    
    % Get spiketimes and shift based on calculated integration time
    
    if exist('Spikes','var')                                 % >>> UMS <<<
        
        spiketimes = unique(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 + spkshift);  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(Clusters(iClu).spikeTimes * 1000 + spkshift)';
        
    end
    
    
    fprintf(' analyzing ch %i clu %i\n',channel,clu)
    
    % Set ymaxval
    yin = find( ( max(yvals, 2+3*max(nanmean(UnitData(iUn).FR_raw_tr,1)) ) - yvals )==0 );
    
    
    %% Get MPH data
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes,RateStream);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
    %%  3) PLOT DATA
    
    T = 250;
    
    for this_rate = AMrates
        
        Period = 1000/this_rate;
        irate     = find(AMrates==this_rate);
        
        MPH_rate  = MPH(MPH.AMrate==this_rate,:);
        
        xvalues  = MPH_rate.Prev100msFR;
%         yvalues  = cellfun(@length, MPH_rate.x) ./ cellfun(@(x) size(x,1), MPH_rate.raster) / Period * 1000;
        yvalues  = cellfun(@(x) mean(x,1), cellfun(@(x) mean(x,2), MPH_rate.FRsmooth, 'UniformOutput',false));
        
        if numel(yvalues)<3
            continue
        end
        [r,p]=corr(xvalues,yvalues);
        
        if p<0.05
%             figure;
%             plot(xvalues,yvalues,'.k')
%             xlim([0 30])
%             ylim([0 15])
%             keyboard
        end
        
        pvalues = [pvalues p];
        
    end  %this_rate
    
end %iUn


figure; hist(pvalues,0:0.01:1)
sum(pvalues<0.05)/length(pvalues)
keyboard

hfcc = figure;
plot([0 1]',[0 1]','-k')
hold on
for irate = 1:3
    plot(Data(iUn,irate).Res_L1o.dprime(:,2),Data(iUn,irate).Res_t10.dprime(:,2),'ok','LineWidth',2)
end
axis square

print_eps_kp(hfcc,fullfile(fn.figs,'MPHclass','CompareClassifiers'))


end %function




