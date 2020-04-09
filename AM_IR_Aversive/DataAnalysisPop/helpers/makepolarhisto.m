function HistOut = makepolarhisto(DegVec,ic,PhaseData)
% HistOut = makepolarhisto(DegVec,ic,PhaseData)
%  
%  Like the matlab function histcounts, but circular.
% 
% KP, 2020-04

HistOut = zeros(numel(DegVec),5);

for irad = 1:numel(DegVec)
    if irad==1 || irad==numel(DegVec)
        flags =  PhaseData(ic,:)>DegVec(end-1) | PhaseData(ic,:)<=DegVec(2);
    else
        flags =  PhaseData(ic,:)>DegVec(irad-1) & PhaseData(ic,:)<=DegVec(irad+1);
    end
    
    HistOut(irad,:) = sum(flags,1);
    
end
end