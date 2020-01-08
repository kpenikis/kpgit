function [CTTS,theseCells,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    whichCells, nTrialMat, UnitData, theseStim, iRS, iNS, minTrs, convwin, AnWin )
% Gets only those values of CTTS to classify.
% KP, 2019-12 



switch whichCells
    case {'all' 'RSNS' 'each'}                                         % All Cells
        theseCells = find(all(nTrialMat(:,theseStim)>=minTrs,2));
        
        if strcmp(whichCells,'RSNS')
            keyboard
            iiRS = intersect(theseCells,iRS);
            iiNS = intersect(theseCells,iNS);
            theseCells = sort([iiRS; iiNS]);
            keyboard
            [~,RS_CTTS] = intersect(theseCells,iiRS);
            [~,NS_CTTS] = intersect(theseCells,iiNS);
        end
        
%     case {'Mar28-AM' 'Mar30-AM' 'Apr02-AM' 'Apr11-AM' 'Jan17-AM' 'Oct26-AM' 'Mar26-AM' 'Jan25-AM' 'Jan21-AM'}
    case unique({UnitData.Session})
        theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,theseStim)>=minTrs,2) );
        
    case 'select'                                  % One or subset of cells
        keyboard
        theseCells = Un_Index;
        
    case 'RS'                                                          % RS
        theseCells = iRS(all(nTrialMat(iRS,theseStim)>=minTrs,2));
        
    case 'NS'                                                          % NS
        theseCells = iNS(all(nTrialMat(iNS,theseStim)>=minTrs,2));
        
%     case {'Best10' 'Worst10'}
%         if nargin<11 || ~exist('theseCells','var')
%             keyboard
%         end
        
        
end


% Filter the data matrix
CTTS = Cell_Time_Trial_Stim(theseCells,AnWin,:,theseStim);


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

