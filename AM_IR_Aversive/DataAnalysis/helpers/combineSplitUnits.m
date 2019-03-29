function combineSplitUnits( UnitData, UnitInfo, SAVENAME )
% Run after AssessUnits. Make sure to manually add unit strings to
% SplitUnits mat file first. Skips units that have already been merged.
% 
%  KP, 2018-06
%


global fn

% Load SplitUnits struct
load([fn.processed '/SplitUnits.mat'])

% Load Units files
if nargin<2
    q = load(fullfile(fn.processed,'Units'));
    UnitData = q.UnitData;
    UnitInfo = q.UnitInfo;
    clear q
end

% Go through each unit and check if it should be merged

for un1_idx = 1:numel(UnitData)
    
    un1_string = [UnitInfo(un1_idx,:).Subject{:}(1:4) '_' UnitInfo(un1_idx,:).Session{:} '_' num2str(UnitInfo(un1_idx,:).Channel) '_' num2str(UnitInfo(un1_idx,:).Clu) ];
    
    % If this string matches one from LUT, proceed
    
    if any(cellfun(@(x) strcmp(un1_string,x), {SplitUnits.un1})) && ~strncmp(UnitInfo.RespType{un1_idx},'merged',6)
        
        un2_string = SplitUnits(cellfun(@(x) strcmp(un1_string,x), {SplitUnits.un1})).un2;
        un2_cell = strsplit(un2_string,'_');
        
        un2_idx = find(strcmp(UnitInfo.Session,un2_cell{2}) & UnitInfo.Channel==str2num(un2_cell{3}) & UnitInfo.Clu==str2num(un2_cell{4}));
        
        if (UnitData(un1_idx).spl ~= UnitData(un2_idx).spl) || (UnitData(un1_idx).lpn ~= UnitData(un2_idx).lpn)
            keyboard
        end
        
        % Add an entry to UnitInto representing joint datapoint
        add_row = {UnitData(un1_idx).Subject ...
            [UnitData(un1_idx).Session UnitData(un2_idx).Session]...
            UnitData(un1_idx).Channel...
            UnitData(un1_idx).Clu * UnitData(un2_idx).Clu...
            UnitInfo.RespType{un1_idx}};
        
        UnitInfo = [UnitInfo; add_row];
        
        UnitInfo.RespType{un1_idx} = 'merged1';
        UnitInfo.RespType{un2_idx} = 'merged2';
        
        
        % Add an entry to UnitData representing joint datapoint
        
        N = numel(UnitData)+1;
        UnitData(N).Subject     = UnitData(un1_idx).Subject;
        UnitData(N).Session     = [UnitData(un1_idx).Session '_' UnitData(un2_idx).Session];
        UnitData(N).Channel     = [UnitData(un1_idx).Channel     UnitData(un2_idx).Channel];
        UnitData(N).Clu         = [UnitData(un1_idx).Clu         UnitData(un2_idx).Clu];
%         UnitData(N).unType      = UnitData(un1_idx).unType;
        UnitData(N).spl         = UnitData(un1_idx).spl;
        UnitData(N).lpn         = UnitData(un1_idx).lpn;
        UnitData(N).BaseFR      = mean([UnitData(un1_idx).BaseFR UnitData(un2_idx).BaseFR]);
        UnitData(N).FR_raw_tr   = [UnitData(un1_idx).FR_raw_tr; UnitData(un2_idx).FR_raw_tr];
        UnitData(N).kw_p        = mean([UnitData(un1_idx).kw_p UnitData(un2_idx).kw_p]);
        UnitData(N).wx_p        = mean([UnitData(un1_idx).wx_p UnitData(un2_idx).wx_p]);
        UnitData(N).VSdata_spk  = reshape( mean([ reshape(UnitData(un1_idx).VSdata_spk,numel(UnitData(un1_idx).VSdata_spk),1)...
            reshape(UnitData(un2_idx).VSdata_spk,numel(UnitData(un2_idx).VSdata_spk),1) ],2) ,3,8);
        UnitData(N).VSdata_gap  = reshape( mean([ reshape(UnitData(un1_idx).VSdata_gap,numel(UnitData(un1_idx).VSdata_gap),1)...
            reshape(UnitData(un2_idx).VSdata_gap,numel(UnitData(un2_idx).VSdata_gap),1) ],2) ,3,8);
        UnitData(N).Phase_spk   = mean([UnitData(un1_idx).Phase_spk; UnitData(un2_idx).Phase_spk],1);
        UnitData(N).Phase_gap   = mean([UnitData(un1_idx).Phase_gap; UnitData(un2_idx).Phase_gap],1);
        UnitData(N).IntTime_spk = mean([UnitData(un1_idx).IntTime_spk UnitData(un2_idx).IntTime_spk]);
        UnitData(N).IntTime_gap = mean([UnitData(un1_idx).IntTime_gap UnitData(un2_idx).IntTime_gap]);
        UnitData(N).ntr         = [UnitData(un1_idx).ntr; UnitData(un2_idx).ntr];
        UnitData(N).FR_nrm      = mean([UnitData(un1_idx).FR_nrm; UnitData(un2_idx).FR_nrm],1,'omitnan');
        UnitData(N).dp_mat      = mean([UnitData(un1_idx).dp_mat(:,2) UnitData(un2_idx).dp_mat(:,2)],2,'omitnan');
        
    end
end


% Save Unit files again
save(fullfile(fn.processed,SAVENAME),'UnitInfo','UnitData','-v7.3');


end
