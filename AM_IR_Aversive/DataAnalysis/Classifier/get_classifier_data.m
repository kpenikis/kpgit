function Res = get_classifier_data(raster,TemplateSize,METRIC,binsize)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each

global Iterations trMin

Res = struct;

rng('shuffle')

% Set binsize depending on decoder type
% switch METRIC
%     case 'FR'
%         bin = min(cellfun(@(x) size(x,2), {raster.raster}));
%     case 'SpV'
%         bin = binsize;
% end

%% Get NoGo (periodic) stimulus info 

NGidx = 1;
NOGO  = raster(NGidx);
rsNOGO  = NOGO.raster;  % Filter to start when sound begins


%% Step through the GO (non-periodic) stimuli to compare to Periodic

Gidxs  = 1+find([raster(2:end).nTrs]>=trMin*2);  %must have more than N trials

% Preallocate
dprime = nan( numel(Gidxs), 1 );
pHit   = nan( numel(Gidxs), 1 );
pFA    = nan( numel(Gidxs), 1 );
nYes   = nan( 1+numel(Gidxs), 1 );
nTrs   = nan( 1+numel(Gidxs), 1 );

for ii = 1:numel(Gidxs)
    
    % Get all identical go trials of this type
    GO     = raster(Gidxs(ii));
    rsGO   = GO.raster;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Run classifier
    [dprime(ii),pHit(ii),pFA(ii)] = run_classifier_Template(rsGO,rsNOGO,Iterations,TemplateSize);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Get data for PY matrix for this GO stim
    [nYes(ii+1),nTrs(ii+1)] = calculate_response_data_phys(0,pHit(ii),pFA(ii),size(rsGO,1),size(rsNOGO,1));
    
end %ii

% Now get data for PY matrix for NOGO stim
[nYes(1),nTrs(1)] = calculate_response_data_phys(1,nan,pFA,nan,size(rsNOGO,1));


% Corrrect inf dprime vals
for inf_idx = find(dprime==Inf)'
    if pHit(inf_idx)==1
        pHit(inf_idx) = 0.999;
    end
    if pFA(inf_idx)==0
        pFA(inf_idx) = 0.001;
    elseif pFA(inf_idx)==1
        pFA(inf_idx) = 0.999;
    end
    corrected_dp = calculate_dprime(pHit(inf_idx),pFA(inf_idx));
    dprime(inf_idx) = corrected_dp;
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Put data into output struct
Res.PYdata   = [ [1 Gidxs]'  nYes  nTrs ];
Res.dprime   = [ Gidxs' dprime];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end