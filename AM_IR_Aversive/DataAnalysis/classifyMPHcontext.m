function Data = classifyMPHcontext(USE_MEASURE)
%
%  classifyMPHcontext
%
%
%
%  KP, 2019-04
%


global fn AMrates rateVec_AC rateVec_DB trMin Iterations

%!!!!!!!!!!!!!!!!!!!
if nargin<1
    USE_MEASURE = 'FR'; 'spikes';
end
%!!!!!!!!!!!!!!!!!!!
RERUN       =   0;
%!!!!!!!!!!!!!!!!!!!
trMin       =  10;
%!!!!!!!!!!!!!!!!!!!
Iterations  =  10;
%!!!!!!!!!!!!!!!!!!!
tVarBin     =  31;
%!!!!!!!!!!!!!!!!!!!
SHUFFLE_SPIKES = 0;
%!!!!!!!!!!!!!!!!!!!


%% Load Unit data files

fn = set_paths_directories;

q = load(fullfile(fn.processed,'Units'));
UnitData = q.UnitData;
UnitInfo = q.UnitInfo;
clear q
%-------
spkshift = mean([UnitData([UnitData.IntTime_spk]>0).IntTime_spk]);
%-------

% Load IR stimulus rate vectors
q = load(fullfile(fn.stim,'rateVec_AC'));
rateVec_AC = q.buffer;
q = load(fullfile(fn.stim,'rateVec_DB'));
rateVec_DB = q.buffer;

AMrates = [2 4 8 16 32];


%% Preallocate

savedir = fullfile(fn.processed,'MPHclassifier');
if ~exist(savedir,'dir')
    mkdir(savedir)
end

if SHUFFLE_SPIKES
    savename = 'ClassData_shuff';
else
    savename = 'ClassData';
end

if RERUN
    Data = struct;
    Un1  = 1;
else
    q = load(fullfile(savedir,savename));
    Data = q.Data;
    clear q
    Un1  = size(Data,1)+1;
end

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

% Make figures
% for irate = 1:numel(AMrates)
%     hf(irate) = figure; hold on
%     plot([AMrates(irate) AMrates(irate)],[0 3],'Color',colors(irate+1,:),'LineWidth',2)
%     set(gca,'xscale','log','xtick',2.^[1:5])
%     xlim(2.^[0.8 5.2])
%     ylim([0 2.5])
%     xlabel('avg AM rate in previous 500 ms')
%     ylabel('d prime')
%     title([num2str(AMrates(irate)) ' Hz period'])
% end



%% Step through Units

for iUn = Un1:numel(UnitData)
        
    %%% skips merged units for now
    if numel(UnitInfo(iUn,:).Session{:})==4  %strncmp(UnitInfo.RespType{iUn},'merged',6)
        continue
    end
    
%     if UnitData(iUn).kw_p>0.05
%         continue
%     end
    
    subject     = UnitData(iUn).Subject;
    session     = UnitData(iUn).Session;
    channel     = UnitData(iUn).Channel(1);
    clu         = UnitData(iUn).Clu(1);
    subjcol     = [1 1 1];
    
    % Get sound parameters
    dBSPL       = UnitData(iUn).spl;
    LP          = UnitData(iUn).lpn;
    
    
    % Load data files
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) )) || iUn==1 || ~exist('TrialData','var')
        fprintf('Loading %s sess %s...\n',subject,session)
        clear TrialData Info RateStream SoundStream SpoutStream
        filename = sprintf( '%s_sess-%s_Info'     ,subject,session); load(fullfile(fn.processed,subject,filename));
        filename = sprintf( '%s_sess-%s_TrialData',subject,session); load(fullfile(fn.processed,subject,filename));
    end
    if (iUn>1 && ~( strcmp(subject,UnitData(iUn-1).Subject) && strcmp(session,UnitData(iUn-1).Session) && channel==UnitData(iUn-1).Channel ) )  || iUn==1 || ~exist('TrialData','var')
        clear Clusters Spikes
        filename = sprintf( '%s_sess-%s_Spikes'   ,subject,session); load(fullfile(fn.processed,subject,filename));
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
    
    
    %%
    if SHUFFLE_SPIKES
        spiketimes = sort(randi( length(SpoutStream), size(spiketimes) ));
    end
    
    
    %% Get MPH data
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    MPH = makeMPHtable(TrialData,Info.artifact(channel).trials',dBSPL,LP,spiketimes,RateStream);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    %%  3) PLOT DATA
    
    try
        
    T = 250;
    
    for this_rate = AMrates
        
        fprintf(' %iHz',this_rate)
        Period = 1000/this_rate;
        irate     = find(AMrates==this_rate);
        
        MPH_rate  = MPH(MPH.AMrate==this_rate,:);
        unique([MPH_rate.PrevAMrt500]);
        
        %~~~~~~~~~~~~~~~~  PERIODIC  ~~~~~~~~~~~~~~~~
        iPdc = find( (MPH_rate.ThisStimID==irate+1 & MPH_rate.Starttime>=T)...
            | (MPH_rate.ThisStimID==irate+1 & MPH_rate.PrevStimID==irate+1) );
        
        if isempty(iPdc)
            fprintf('   ...skipping bc no valid Pdc MPH\n')
            continue
        end
        
        nTrs = cellfun(@(x) size(x,1),MPH_rate(iPdc,:).raster);
        thisMPH = struct();
        np=1;
        
        thisMPH(1).Context      = 'Pdc'; 
        thisMPH(1).indices      = iPdc;
        switch USE_MEASURE
            case 'FR'
                thisMPH(1).raster  = vertcat(MPH_rate(iPdc,:).FRsmooth{:});
            case 'spikes'
                thisMPH(1).raster  = vertcat(MPH_rate(iPdc,:).raster{:});
        end
        thisMPH(1).Prev500msFR  = cell2mat(MPH_rate(iPdc,:).Prev500msFR')'; 
        thisMPH(1).Prev100msFR  = cell2mat(MPH_rate(iPdc,:).Prev100msFR')'; 
        thisMPH(1).PrevAMrt500  = sum([MPH_rate(iPdc,:).PrevAMrt500] .* nTrs) / sum(nTrs);
        thisMPH(1).PrevAMrt100  = sum([MPH_rate(iPdc,:).PrevAMrt100] .* nTrs) / sum(nTrs);
        thisMPH(1).nTrs         = sum(nTrs);
        
        
        %~~~~~~~~~~~~~~~~  IRREGULAR  ~~~~~~~~~~~~~~~~
        
        theseIRs = unique(MPH_rate.ThisStimID(MPH_rate.ThisStimID>6))';
        
        for thisIR = theseIRs
            
            PPds = unique(MPH_rate.PrevPd( MPH_rate.Starttime>=T & MPH_rate.ThisStimID==thisIR ));
            
            for thisPP = PPds'
                
                iIR = find(MPH_rate.PrevPd==thisPP & MPH_rate.Starttime>T & MPH_rate.ThisStimID==thisIR)';
                if isempty(iIR), continue, end
                
                nTrs = cellfun(@(x) size(x,1),MPH_rate(iIR,:).raster);
                
                np = np+1;
                thisMPH(np).Context      = 'IR';
                thisMPH(np).indices      = iIR;
                
                switch USE_MEASURE
                    case 'FR'
                        thisMPH(np).raster  = vertcat(MPH_rate(iIR,:).FRsmooth{:});
                    case 'spikes'
                        thisMPH(np).raster  = vertcat(MPH_rate(iIR,:).raster{:});
                end
                
                thisMPH(np).Prev500msFR  = cell2mat(MPH_rate(iIR,:).Prev500msFR')';
                thisMPH(np).Prev100msFR  = cell2mat(MPH_rate(iIR,:).Prev100msFR')';
                thisMPH(np).PrevAMrt500  = sum([MPH_rate(iIR,:).PrevAMrt500] .* nTrs) / sum(nTrs);
                thisMPH(np).PrevAMrt100  = sum([MPH_rate(iIR,:).PrevAMrt100] .* nTrs) / sum(nTrs);
                thisMPH(np).nTrs         = sum(nTrs);
                
            end %PPds
        end %thisIR
        
        
        %~~~~~~~~~~~~~~~~  Post-Transition  ~~~~~~~~~~~~~~~~
        
        iTrans = find(MPH_rate.ThisStimID==irate+1 & MPH_rate.Starttime<T/2 & MPH_rate.PrevStimID~=irate+1 )';
        PStIDs = [MPH_rate(iTrans,:).PrevStimID]';
        
        for ips = unique(PStIDs)
            
            np=np+1;
            iTr = iTrans(PStIDs==ips);
            
            nTrs                     = cellfun(@(x) size(x,1),MPH_rate(iTr,:).raster);
            
            thisMPH(np).Context      = 'Trans';
            thisMPH(np).indices      = iTr;
            switch USE_MEASURE
                case 'FR'
                    thisMPH(np).raster  = vertcat(MPH_rate(iTr,:).FRsmooth{:});
                case 'spikes'
                    thisMPH(np).raster  = vertcat(MPH_rate(iTr,:).raster{:});
            end
            thisMPH(np).Prev500msFR  = cell2mat(MPH_rate(iTr,:).Prev500msFR')';
            thisMPH(np).Prev100msFR  = cell2mat(MPH_rate(iTr,:).Prev100msFR')';
            thisMPH(np).PrevAMrt500  = sum([MPH_rate(iTr,:).PrevAMrt500] .* nTrs) / sum(nTrs);
            thisMPH(np).PrevAMrt100  = sum([MPH_rate(iTr,:).PrevAMrt100] .* nTrs) / sum(nTrs);
            thisMPH(np).nTrs         = sum(nTrs);
            
        end %PStIDs
        
        
        % Run classifiers and save data for later
        Data(iUn,irate).AMrate   = this_rate;
        Data(iUn,irate).data     = thisMPH;
        fprintf('   running classifier: leave one tr out\n')
        Data(iUn,irate).Res_L1o  = get_classifier_data( thisMPH, -1 );
%         fprintf('   running classifier: 10 trial templates\n')
        Data(iUn,irate).Res_t10  = []; %get_classifier_data( thisMPH, 10 );
        
        fprintf(' \n')
        
        % Add to plot
%         these_idx       = Data(iUn,irate).Res_t10.dprime(:,1);
%         x_plot          = [Data(iUn,irate).data(these_idx).PrevAMrt500];
%         [x_plot,x_idx]  = sort(x_plot);
%         
%         figure(hf(irate));
%         plot(x_plot, Data(iUn,irate).Res_t10.dprime(x_idx,2)',...
%             '-ok')
%
%         if any([Data(iUn,irate).Res_t10.dprime(:,2)]>1)
%             aaaa=234;
%         end
        
    end  %this_rate
    
    
    % Save Data to tmp folder
    save(fullfile(savedir,'tmp',savename),'Data','-v7.3')
    fprintf('\n')
    
    catch
        keyboard
    end
    
end %iUn

% Save final Data struct
save(fullfile(savedir,savename),'Data','-v7.3')

return


hfcc = figure;
plot([0 1]',[0 1]','-k')
hold on
for irate = 1:5
    plot(Data(iUn,irate).Res_L1o.dprime(:,2),Data(iUn,irate).Res_t10.dprime(:,2),'ok','LineWidth',2)
end
axis square

print_eps_kp(hfcc,fullfile(fn.figs,'MPHclass','CompareClassifiers'))


end %function




