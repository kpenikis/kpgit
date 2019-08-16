function Data = ClassifyMPHcontext_matched(USE_MEASURE)
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
RERUN       =   1;
%!!!!!!!!!!!!!!!!!!!
trMin       =  10;
%!!!!!!!!!!!!!!!!!!!
Iterations  =  1000;
%!!!!!!!!!!!!!!!!!!!
tVarBin     =  31;
%!!!!!!!!!!!!!!!!!!!
SHUFFLE_SPIKES = 1;
%!!!!!!!!!!!!!!!!!!!
T         = 250;
alfa      = 0.05;
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
    q = load(fullfile(savedir,'tmp',savename));
    Data = q.Data;
    clear q
    Un1  = size(Data,1)+1;
end


%% Step through Units

for iUn = Un1:numel(UnitData)
    
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
        clear TrialData Info RateStream SoundStream SpoutStream RateStream
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
        
        spiketimes = unique(Spikes.sorted(channel).spiketimes(Spikes.sorted(channel).assigns==clu') * 1000 - spkshift);  %ms
        
    elseif exist('Clusters','var')                            % >>> KS <<<
        
        iClu = find([Clusters.maxChannel] == channel & [Clusters.clusterID] == clu);
        spiketimes = unique(Clusters(iClu).spikeTimes * 1000 - spkshift)';
        
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
    
    
    
    %%  3) GET MATCHED PERIODS, COLLECT DATA
    
    
    for this_rate = AMrates   %(UnitData(iUn).iBMF_FR)
        
        fprintf(' %iHz',this_rate)
        
        % Set up
        thisMPH   = struct;
        Period    = 1000/this_rate;
        irate     = find(AMrates==this_rate);
        
        MPH_rate  = MPH(MPH.AMrate==this_rate,:);
        
%         % Get Pdc response type                      [ early  late  p_val ]
%         PRT = 'na';
%         if ~isempty(UnitData(iUn).DeltaNspk) && ~isempty(UnitData(iUn).DeltaNspk{this_rate==AMrates})
%             if (UnitData(iUn).DeltaNspk{this_rate==AMrates}(3))<alfa && (UnitData(iUn).DeltaNspk{this_rate==AMrates}(1) > UnitData(iUn).DeltaNspk{this_rate==AMrates}(2))
%                 % adapting
%                 PRT = 'A';
%             elseif (UnitData(iUn).DeltaNspk{this_rate==AMrates}(3))<alfa && (UnitData(iUn).DeltaNspk{this_rate==AMrates}(1) < UnitData(iUn).DeltaNspk{this_rate==AMrates}(2))
%                 % facilitating
%                 PRT = 'F';
%             else
%                 PRT = 'NC';
%             end
%         end
        
        np=0;
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Get instances of this period in IRREGULAR context first
        
        theseIRs = unique(MPH_rate.ThisStimID(MPH_rate.ThisStimID>6))';
        
        for thisIR = theseIRs
            
            for iseq = 1:2
                
                iIR = find( MPH_rate.ThisStimID==thisIR & MPH_rate.ThisStimID==thisIR & MPH_rate.SeqPos==iseq )';
                
                if isempty(MPH_rate(iIR,:)), continue, end
                
                np=np+1;
                
                nTrs = cellfun(@(x) size(x,1),MPH_rate(iIR,:).raster);
                
                thisMPH(np,2).Context      = 'IR';
                thisMPH(np,2).indices      = iIR;
                
                switch USE_MEASURE
                    case 'FR'
                        thisMPH(np,2).raster  = vertcat(MPH_rate(iIR,:).FRsmooth{:});
                    case 'spikes'
                        thisMPH(np,2).raster  = vertcat(MPH_rate(iIR,:).raster{:});
                end
                
                thisMPH(np,2).Prev500msFR  = cell2mat(MPH_rate(iIR,:).Prev500msFR')';
                thisMPH(np,2).Prev100msFR  = cell2mat(MPH_rate(iIR,:).Prev100msFR')';
                thisMPH(np,2).PrevAMrt500  = sum([MPH_rate(iIR,:).PrevAMrt500] .* nTrs) / sum(nTrs);
                thisMPH(np,2).PrevAMrt100  = sum([MPH_rate(iIR,:).PrevAMrt100] .* nTrs) / sum(nTrs);
                thisMPH(np,2).nTrs         = sum(nTrs);
                
            end %iseq
            
        end %IRstim -- IRREGULAR first
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Now get matching PERIODIC context data
        
        % Get Pdc starttimes, but SKIP if too close to trial transition
        PdcMPs = MPH_rate(MPH_rate.ThisStimID<7,:);
        pdcstarttimes = unique(PdcMPs.Starttime)';
        pdcstarttimes(pdcstarttimes<T) = nan;
        
        
        % Now get a matching Pdc MP for each Irr one
        
        for ip = 1:size(thisMPH,1)
            
            % Find matching period by closest starttime
            starttime = mode(MPH_rate(thisMPH(ip,2).indices,:).Starttime);
            
            strtdiffs = abs(starttime - pdcstarttimes);
            [STdiff,ipd] = min(strtdiffs);
            if STdiff>1000
                aaa=234;
            end
            
            iPdc = find((MPH_rate.ThisStimID<7 & MPH_rate.Starttime==pdcstarttimes(ipd)))';
            
            if isempty(iPdc)
                thisMPH(ip,1).Context      = 'Pdc';
                thisMPH(ip,1).indices      = iPdc;
                thisMPH(ip,1).raster       = [];
                thisMPH(ip,1).Prev500msFR  = [];
                thisMPH(ip,1).Prev100msFR  = [];
                thisMPH(ip,1).PrevAMrt500  = [];
                thisMPH(ip,1).PrevAMrt100  = [];
                thisMPH(ip,1).nTrs         = 0;
                
            else
                nTrs = cellfun(@(x) size(x,1),MPH_rate(iPdc,:).raster);
                
                thisMPH(ip,1).Context      = 'Pdc';
                thisMPH(ip,1).indices      = iPdc;
                switch USE_MEASURE
                    case 'FR'
                        thisMPH(ip,1).raster  = vertcat(MPH_rate(iPdc,:).FRsmooth{:});
                    case 'spikes'
                        thisMPH(ip,1).raster  = vertcat(MPH_rate(iPdc,:).raster{:});
                end
                thisMPH(ip,1).Prev500msFR  = cell2mat(MPH_rate(iPdc,:).Prev500msFR')';
                thisMPH(ip,1).Prev100msFR  = cell2mat(MPH_rate(iPdc,:).Prev100msFR')';
                thisMPH(ip,1).PrevAMrt500  = sum([MPH_rate(iPdc,:).PrevAMrt500] .* nTrs) / sum(nTrs);
                thisMPH(ip,1).PrevAMrt100  = sum([MPH_rate(iPdc,:).PrevAMrt100] .* nTrs) / sum(nTrs);
                thisMPH(ip,1).nTrs         = sum(nTrs);
            end
        end %ip
        
        
        %% Run classifiers and save data for later
        
        Data(iUn,irate).AMrate   = this_rate;
        Data(iUn,irate).data     = thisMPH;
        fprintf('   running classifier: leave one tr out\n')
        Data(iUn,irate).Res_L1o  = get_classifier_data( thisMPH, -1 );
        %         fprintf('   running classifier: 10 trial templates\n')
        Data(iUn,irate).Res_t10  = []; %get_classifier_data( thisMPH, 10 );
        
        fprintf(' \n')
        
        
    end  %this_rate
    
    
    %% Save Data to tmp folder
    
    save(fullfile(savedir,'tmp',savename),'Data','-v7.3')
    fprintf('\n')
    
    
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




