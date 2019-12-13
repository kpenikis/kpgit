function [CTTS,nUns,Dur,nStim] = filterDataMatrix( Cell_Time_Trial_Stim, ...
    whichIrr, whichCells, nTrialMat, UnitData, ACstim, DBstim, iRS, iNS, minTrs, convwin )

switch whichCells
    case {'all' 'RSNS'}                                         % All Cells
        switch whichIrr
            case 'AC'
                theseCells = find(all(nTrialMat(:,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = find(all(nTrialMat(:,DBstim)>=minTrs,2));
        end
        
        if strcmp(whichCells,'RSNS')
            keyboard
            iiRS = intersect(theseCells,iRS);
            iiNS = intersect(theseCells,iNS);
            theseCells = sort([iiRS; iiNS]);
            keyboard
            [~,RS_CTTS] = intersect(theseCells,iiRS);
            [~,NS_CTTS] = intersect(theseCells,iiNS);
        end
        
        
    case {'Mar28-AM' 'Mar30-AM' 'Apr02-AM' 'Apr11-AM'}            % Session
        
        switch whichIrr
            case 'AC'
                theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,ACstim)>=minTrs,2) );
            case 'DB'
                theseCells = find( strcmp({UnitData.Session}',whichCells) & all(nTrialMat(:,DBstim)>=minTrs,2) );
        end
        
        
    case 'select'                                  % One or subset of cells
        keyboard
        theseCells = Un_Index;
        
        
    case 'RS'                                                          % RS
        switch whichIrr
            case 'AC'
                theseCells = iRS(all(nTrialMat(iRS,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = iRS(all(nTrialMat(iRS,DBstim)>=minTrs,2));
        end
        
        
    case 'NS'                                                          % NS
        switch whichIrr
            case 'AC'
                theseCells = iNS(all(nTrialMat(iNS,ACstim)>=minTrs,2));
            case 'DB'
                theseCells = iNS(all(nTrialMat(iNS,DBstim)>=minTrs,2));
        end
end


switch whichIrr
    case 'AC'
        CTTS = Cell_Time_Trial_Stim(theseCells,:,:,ACstim);
        nTrialMat = nTrialMat(theseCells,ACstim);
    case 'DB'
        CTTS = Cell_Time_Trial_Stim(theseCells,:,:,DBstim);
        nTrialMat = nTrialMat(theseCells,DBstim);
end

% Repmat, to artificially increase dimensionality of training data
CTTS = repmat(CTTS,[50 1 1 1]);


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

