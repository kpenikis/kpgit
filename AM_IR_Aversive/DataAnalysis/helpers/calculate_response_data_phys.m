
function [pYes,nTrs] = calculate_response_data_phys(ngL,pHit,pFA,nGO,nNG)
%  Calculate response data based on classifier output.
%
%  Inputs
%    ngL:  logical, True if called for nogo
%    pHit: hit rate from classifier output for this go stim
%    pFA:  FA rate from classifier output for all go stimuli
%    nGO:  number of trials for this go stim
%    nNG:  number of trials for the nogo stim
% 
%  KP, 2017-02
% 

%NOGO stimulus
if ngL==1 
    
    % Get average FA rate
    pFA = mean(pFA(2:end));
    
    %Correct floor
    if pFA <0.05
        pFA = 0.05;
    end
    %Correct ceiling
    if pFA >0.95
        pFA = 0.95;
    end
    
    pYes = pFA*nNG;
    nTrs = nNG;

%GO stimulus
else
    
    %Correct floor
    if pHit <0.05
        pHit = 0.05;
    end
    %Correct ceiling
    if pHit >0.95
        pHit = 0.95;
    end
    
    pYes = pHit*nGO;
    nTrs = nGO;
    
end


end


