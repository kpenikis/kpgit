function [i_sorted,sortdata] = sort_thLat_BegEnd(data)
% called by PopMPH 
%

global Boundaries


sortdata = [];

% do it the slow but easy way
for iUn=1:size(data,1)
    
    clear uTh uLat uAbv
    
    if all(isnan(data(iUn,:)))
        sortdata = [sortdata; nan nan nan nan ];
    else
        
        %find highest thresh exceeded
        uTh = Boundaries(find(max(data(iUn,:),[],2) > Boundaries,1,'last'));
        %find time exceeded ith
        uLatB = find(data(iUn,:) > uTh,1,'first');
        %find time crossed below ith
        uLatE = size(data,2) - find(data(iUn,:) > uTh,1,'last');
        %find ms above
        uAbv = sum(data(iUn,:) > uTh);
        
        % Save unit data
        sortdata = [sortdata; uTh round(uLatB/5)*5 round(uLatE/5)*5 uAbv ];
    end
end
[~,i_sorted] = sortrows( sortdata, [1 -2 3 4]);


end

