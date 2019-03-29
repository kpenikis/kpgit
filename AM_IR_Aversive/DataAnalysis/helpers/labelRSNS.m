function UnitInfo = labelRSNS(UnitInfo)
% Adds field SpkShape to UnitInfo to label each unit as narrow spiking
% (NS) or regular spiking (RS). Current rule is that NS units have a width
% at half max <=0.23 ms AND time from trough to peak <0.5 ms. 
%
% KP, 2019-03
% 

% Get indices of narrow spikes and regular spiking units
SpkShape = cell(size(UnitInfo,1),1);
SpkShape([UnitInfo.WidthHalfMax]<=0.23 & [UnitInfo.TroughPeak]<0.5) = deal({'NS'});
SpkShape([UnitInfo.WidthHalfMax]>0.23 & [UnitInfo.TroughPeak]>=0.5) = deal({'RS'});

UnitInfo.SpkShape = SpkShape;

end