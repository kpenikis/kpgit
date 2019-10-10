function [i_sorted,sortdata] = sort_thLat(data)
% called by PopMPH 
%

global Boundaries


sortdata = [];

% do it the slow but easy way
for iUn=1:size(data,1)
    
    clear uTh uLat uAbv
    
    if all(isnan(data(iUn,:)))
        sortdata = [sortdata; nan nan nan ];
    else
        
        %find highest thresh exceeded
        uTh = Boundaries(find(max(data(iUn,:),[],2) > Boundaries,1,'last'));
        %find time exceeded it
        uLat = find(data(iUn,:) > uTh,1,'first');
        %find ms above
        uAbv = sum(data(iUn,:) > uTh);
        
        % Save unit data
        sortdata = [sortdata; uTh uLat uAbv ];
    end
end
[~,i_sorted] = sortrows( sortdata, [1 -2 -3]);


end

