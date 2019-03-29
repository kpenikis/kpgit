function FRmat = format_FRmat( NOGOdata, GOdata )
% called by AssessResponses
% FRmat = [ stim_id meanFR stdFR semFR ]
% KP, 2018-03


FRmat = nan(size(GOdata,3)+1,4);

% Do NOGO data first 
meanFR = mean(sum(NOGOdata,2)/size(NOGOdata,2)*1000);
stdFR  = std(sum(NOGOdata,2)/size(NOGOdata,2)*1000);
semFR  = stdFR/sqrt(size(NOGOdata,1)-1);

FRmat(1,:) = [0 meanFR stdFR semFR];


% Now get data from GO stimuli
for is = 1:size(GOdata,3)
    
    data = GOdata(:,:,is);
    
    meanFR = mean(sum(data,2)/size(data,2)*1000);
    stdFR  = std(sum(data,2)/size(data,2)*1000);
    semFR  = stdFR/sqrt(size(data,1)-1);
        
    FRmat(is+1,:) = [is meanFR stdFR semFR];
    
%     if any(isnan(FRmat(is+1,:)))
%         keyboard
%     end
end


end