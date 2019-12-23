function [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    whichCells, nTrialMat, UnitData, theseStim, iRS, iNS, minTrs, convwin, AnWin )

switch whichCells
    case {'all' 'RSNS'}                                         % All Cells
        theseCells = find(all(nTrialMat(:,theseStim)>=minTrs,2));
%         switch whichIrr
%             case 'AC'
%                 theseCells = find(all(nTrialMat(:,ACstim)>=minTrs,2));
%             case 'DB'
%                 theseCells = find(all(nTrialMat(:,DBstim)>=minTrs,2));
%             case 'Speech'
%                 theseCells = find(all(nTrialMat>=minTrs,2));
%         end
        
        if strcmp(whichCells,'RSNS')
            keyboard
            iiRS = intersect(theseCells,iRS);
            iiNS = intersect(theseCells,iNS);
            theseCells = sort([iiRS; iiNS]);
            keyboard
            [~,RS_CTTS] = intersect(theseCells,iiRS);
            [~,NS_CTTS] = intersect(theseCells,iiNS);
        end
        
        
    case {'Mar28-AM' 'Mar30-AM' 'Apr02-AM' 'Apr11-AM' 'Jan17-AM' 'Oct26-AM' 'Mar26-AM' 'Jan25-AM' 'Jan21-AM'}            % Session
        
        theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,theseStim)>=minTrs,2) );
%         switch whichIrr
%             case 'AC'
%                 theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,ACstim)>=minTrs,2) );
%             case 'DB'
%                 theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,DBstim)>=minTrs,2) );
%         end
        
        
    case 'select'                                  % One or subset of cells
        keyboard
        theseCells = Un_Index;
        
        
    case 'RS'                                                          % RS
        theseCells = iRS(all(nTrialMat(iRS,theseStim)>=minTrs,2));
%         switch whichIrr
%             case 'AC'
%                 theseCells = iRS(all(nTrialMat(iRS,ACstim)>=minTrs,2));
%             case 'DB'
%                 theseCells = iRS(all(nTrialMat(iRS,DBstim)>=minTrs,2));
%             case 'Speech'
%                 theseCells = iRS(all(nTrialMat(iRS,:)>=minTrs,2));
%         end
        
        
    case 'NS'                                                          % NS
        theseCells = iRS(all(nTrialMat(iRS,theseStim)>=minTrs,2));
%         switch whichIrr
%             case 'AC'
%                 theseCells = iNS(all(nTrialMat(iNS,ACstim)>=minTrs,2));
%             case 'DB'
%                 theseCells = iNS(all(nTrialMat(iNS,DBstim)>=minTrs,2));
%             case 'Speech'
%                 theseCells = iNS(all(nTrialMat(iNS,:)>=minTrs,2));
%         end
end


% Filter the data matrix
CTTS = Cell_Time_Trial_Stim(theseCells,AnWin,:,theseStim);

% switch whichIrr
%     case 'AC'
%         CTTS = Cell_Time_Trial_Stim(theseCells,AnWin,:,ACstim);
% %         nTrialMat = nTrialMat(theseCells,ACstim);
%     case 'DB'
%         CTTS = Cell_Time_Trial_Stim(theseCells,AnWin,:,DBstim);
% %         nTrialMat = nTrialMat(theseCells,DBstim);
%     case 'Speech'
%         CTTS = Cell_Time_Trial_Stim(theseCells,AnWin,:,:);
% end

% Repmat, to artificially increase dimensionality of training data
% CTTS = repmat(CTTS,[50 1 1 1]);


% Complete
nUns      = size(CTTS,1);
Dur       = size(CTTS,2);
nStim     = size(CTTS,4);


% Convolve now, so don't have to do it within the bootstrapping loop
for i1 = 1:size(CTTS,1)
    for i3 = 1:size(CTTS,3)
        for i4 = 1:size(CTTS,4)
            conv_data = conv(CTTS(i1,:,i3,i4),convwin);
            CTTS(i1,:,i3,i4) = conv_data(1:Dur);
        end
    end
end

end

