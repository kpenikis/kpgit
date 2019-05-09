function StructOut = getFRhist( StructIn, UnitData, spkshift )

global AMrates

Subject = UnitData.Subject;
Session = UnitData.Session;
Channel = UnitData.Channel;
Clu     = UnitData.Clu;
dBSPL   = UnitData.spl;
LP      = UnitData.lpn;

fn = set_paths_directories(Subject,Session);

filename = sprintf( '%s_sess-%s_Info'     ,Subject,Session); load(fullfile(fn.processed,Subject,filename));
filename = sprintf( '%s_sess-%s_TrialData',Subject,Session); load(fullfile(fn.processed,Subject,filename));
filename = sprintf( '%s_sess-%s_Spikes'   ,Subject,Session); load(fullfile(fn.processed,Subject,filename));


% Get spiketimes and shift based on calculated integration time
if exist('Spikes','var')                                 % >>> UMS <<<
    
    spiketimes = unique(Spikes.sorted(Channel).spiketimes(Spikes.sorted(Channel).assigns==Clu') * 1000 + spkshift);  %ms
    
elseif exist('Clusters','var')                           %  >>> KS <<<
    
    iClu = find([Clusters.maxChannel] == Channel & [Clusters.clusterID] == Clu);
    spiketimes = unique(Clusters(iClu).spikeTimes * 1000 + spkshift)';
    
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
MPH = makeMPHtable(TrialData,Info.artifact(Channel).trials',dBSPL,LP,spiketimes,RateStream);
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|||<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



StructOut = StructIn;
for ir = 1:size(StructIn,2)
    for id = 1:size(StructIn(1,ir).data,2)
        
        % These indices in MPH table
        idx = StructIn(1,ir).data(id).indices;
        idx_rate  = find(MPH.AMrate==AMrates(ir));
        idx = idx_rate(idx);
        if iscolumn(idx)
            idx = idx';
        end
        
        % Add full vector of FR history data
        StructOut(1,ir).data(id).Prev500msFR = cell2mat(MPH(idx,:).Prev500msFR')';
        StructOut(1,ir).data(id).Prev100msFR = cell2mat(MPH(idx,:).Prev100msFR')';
        
        % Also update raster and nTrs to exclude skipped trials
        try
            StructOut(1,ir).data(id).raster = vertcat(cell2mat(MPH(idx,:).raster));
            StructOut(1,ir).data(id).nTrs   = size(StructOut(1,ir).data(id).raster,1);
        catch
            keyboard
        end
        
        if size(cell2mat(MPH(idx,:).Prev500msFR'),2) ~= StructOut(1,ir).data(id).nTrs
            keyboard
        end
                
    end
end



end %function



