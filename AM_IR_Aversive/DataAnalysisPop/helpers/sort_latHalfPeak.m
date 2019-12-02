function [i_sorted,sortdata] = sort_latHalfPeak(data)
% called by PopMPH 
%
%       NEVER FINISHED -- MAY BE DEAD END
%

% global Boundaries


sortdata = [];

% do it the slow but easy way
for iUn=1:size(data,1)
    
    clear uTh uLat uAbv
    
    if all(isnan(data(iUn,:)))
        sortdata = [sortdata; nan nan nan nan];
    else
        
        %define half dyn range
        uTh = (max(data(iUn,:))-min(data(iUn,:)))/2+min(data(iUn,:));
        
        [PKS,LOCS] = findpeaks(data(iUn,:),'MinPeakHeight',uTh,'MinPeakWidth',15);
        
        uLat_1p  = LOCS(1);
        uHt_1p   = PKS(1);
        uLat_pM  = LOCS(PKS==max(PKS));
        uHt_pM   = max(PKS);
        
%         %find highest thresh exceeded
%         uTh = Boundaries(find(max(data(iUn,:),[],2) > Boundaries,1,'last'));
%         %find time exceeded it
%         uLat = find(data(iUn,:) > uTh,1,'first');
%         %find ms above
%         uAbv = sum(data(iUn,:) > uTh);
        
        % Save unit data
        sortdata = [sortdata; uLat_1p uHt_1p uLat_pM uHt_pM ];
    end
end
[~,i_sorted] = sortrows( sortdata, [-1 2]);


end

