function Res = get_classifier_data(Data,TemplateSize)
% output = get_classifier_data(raster,TemplateSize)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each

global Iterations trMin

rng('shuffle')


% Preallocate

dprime = nan(size(Data,1),1);
pHit   = nan(size(Data,1),1);
pFA    = nan(size(Data,1),1);
nYes   = nan(size(Data,1),2);
nTrs   = nan(size(Data,1),2);


for ip = 1:size(Data,1)
    
    % Check that both MPs have enough trials
    if ~all([Data(ip,:).nTrs]>=(trMin*2))
        continue
    end
    
    
    % Get NoGo (Pdc) stimulus data
    NOGO   = Data(ip,1);
    rsNOGO = NOGO.raster;
    
    
    % Get GO (Irr) stimulus data
    GO     = Data(ip,2);
    rsGO   = GO.raster;
    
    
    % Skip if no activity at all
    if sum(sum(rsGO))==0 && sum(sum(rsNOGO))==0
        continue
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Run classifier
    [dprime(ip),pHit(ip),pFA(ip)] = run_classifier_Template(rsGO,rsNOGO,Iterations,TemplateSize);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    % Get data for PY matrix for this GO stim
    [nYes(ip,2),nTrs(ip,2)] = calculate_response_data_phys(0,pHit(ip),pFA(ip),size(rsGO,1),size(rsNOGO,1));
    
    % Now get data for PY matrix for NOGO stim
    [nYes(ip,1),nTrs(ip,1)] = calculate_response_data_phys(1,nan,pFA,nan,size(rsNOGO,1));
    
    
end %ip (each MP pair)


% Corrrect inf dprime vals
for inf_idx = find(dprime==Inf)'
    if pHit(inf_idx)==1
        pHit(inf_idx) = 0.99;
    elseif pHit(inf_idx)==0
        pHit(inf_idx) = 0.01;
    end
    if pFA(inf_idx)==0
        pFA(inf_idx) = 0.01;
    elseif pFA(inf_idx)==1
        pFA(inf_idx) = 0.99;
    end
    corrected_dp = calculate_dprime(pHit(inf_idx),pFA(inf_idx));
    dprime(inf_idx) = abs(corrected_dp);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Put data into output struct
Res = struct;
Res.pHpFA    = [ (1:size(Data,1))'  pHit  pFA ];
Res.dprime   = [ (1:size(Data,1))'  dprime];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end


