function output = get_classifier_data(raster,METRIC,binsize,indVar)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each

%%%%%%%%%%%%%%%%%%
iterations = 1000;
%%%%%%%%%%%%%%%%%%


% Get unique depths and jitter vectors
[depths,~,ndpth] = unique([raster.AMdepth]);
[fIDs,~,nfid] = unique([raster.jitter],'stable');

% Set independent and condition variables
[NGval,indValues,indIdx,indLabels,condIdx,condLabels] = ...
    set_variables(indVar,depths,ndpth,fIDs,nfid);

if numel(indValues)<3  %skip data that can't yield a depth function
    output = 'not enough stimuli for neurometric data';
    return
end


% Get minimum stimulus duration
mdur = min([raster.stimDur]);



%% Get nogo stimulus info

NGidx = indIdx==find(indValues==NGval);
NOGO = raster(NGidx);
if numel(NOGO)>1
    if strcmp(indVar,'jitter') && any([NOGO.AMdepth]==0)
        NOGO([NOGO.AMdepth]==0)=[];
    else
        keyboard
        NOGO = collapse_blocks(NOGO);
    end
elseif numel(NOGO)<1
    keyboard
end

% Get nogo spiketimes, starting when sound begins
nogo_x = NOGO.x(NOGO.x>0);
nogo_y = NOGO.y(NOGO.x>0);


% Set up dprime matrix
dprime = nan( numel(condLabels), numel(indValues) );


%% Step through the go stimuli to compare to nogo

for ic = 1:max(condIdx)
    
    if sum(condIdx==ic)<3, continue, end
    
    pHit = nan( numel(indValues), 1 );
    pFA  = nan( numel(indValues), 1 );
    nYes = nan( numel(indValues), 1 );
    nTrs = nan( numel(indValues), 1 );
    
    for ii = 1:max(indIdx)
        
        % Skip depth of 0 (nogo) and fidXdepth conditions that
        % weren't presented.
        if ii==1 && ( (strcmp(indVar,'jitter')&&indValues(ii)==NGval) || (strcmp(indVar,'depth')&&indValues(ii)==convert_depth_proptodB(NGval)) )
            continue
        elseif ii~=1 && ( (strcmp(indVar,'jitter')&&indValues(ii)==NGval) || (strcmp(indVar,'depth')&&indValues(ii)==convert_depth_proptodB(NGval)) )
            keyboard
        end
        if numel(raster((condIdx==ic)&(indIdx==ii)))<1, continue, end
        
        
        % Get all identical go trials of this type
        GO = collapse_blocks(raster((condIdx==ic)&(indIdx==ii)));
        
        % Get spiketimes for go stimulus, starting when sound begins
        try
            go_x = GO.x(GO.x>0);
            go_y = GO.y(GO.x>0);
        catch
            keyboard
        end
        
        % Clip spiketimes at end of stimulus
        nogo_y = nogo_y(nogo_x<=mdur);
        nogo_x = nogo_x(nogo_x<=mdur);
        go_y   = go_y(go_x<=mdur);
        go_x   = go_x(go_x<=mdur);
        
        % Set binsize depending on decoder type
        switch METRIC
            case 'FR'
                bin = mdur;
            case 'SpV'
                bin = binsize;
        end
        
        % Convert raster format to accomodate distance calculations
        GOs = zeros(max(go_y),mdur);
        rsGO = zeros(max(go_y),floor(mdur/bin));
        for iy = 1:max(go_y)
            GOs(iy,go_x(go_y==iy)) = 1;
            rsGO(iy,:) = sum( reshape( GOs(iy,1:(bin*floor(mdur/bin))) , [bin,floor(mdur/bin)] ) ,1);
        end
        
        NOGOs = zeros(max(nogo_y),mdur);
        rsNOGO = zeros(max(nogo_y),floor(mdur/bin));
        for iy = 1:max(nogo_y)
            NOGOs(iy,nogo_x(nogo_y==iy)) = 1;
            rsNOGO(iy,:) = sum( reshape( NOGOs(iy,1:(bin*floor(mdur/bin))) , [bin,floor(mdur/bin)] ),1);
        end
        
        
        %%
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Run classifier
        [dprime(ic,ii),pHit(ii),pFA(ii)] = run_classifier(rsGO,rsNOGO,iterations);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % Get data for PY matrix for this GO stim
        [nYes(ii,1),nTrs(ii,1)] = calculate_response_data_phys(0,pHit(ii),pFA,size(rsGO,1),size(rsNOGO,1));
        
        
    end %ii
    
    
    % Now get data for PY matrix for NOGO stim
    [nYes(1,1),nTrs(1,1)] = calculate_response_data_phys(1,nan,pFA,nan,size(rsNOGO,1));
    
    
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Put data into output struct
    output.PYdata{ic}   = [ indValues'  nYes  nTrs ];
    output.dprime{ic}   = [ indValues(2:end)' dprime(ic,2:end)'];
    output.stim{ic,1}   = condLabels{ic};
    output.stim{ic,2}   = indLabels';
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
end %ic




end