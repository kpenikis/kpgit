function PhysOut = pp_common_avg_ref_data( PhysIn, useChs )
% function PhysOut = pp_common_avg_ref_data( PhysIn, useChs )
% Step through each channel and subtract the mean of all other channels. 
% KP,2018-12
%

PhysOut = PhysIn;

for ich = useChs
    theseChs = useChs(useChs~=ich);
    PhysOut(ich,:) = PhysIn(ich,:) - mean(PhysIn(theseChs,:));
end

end