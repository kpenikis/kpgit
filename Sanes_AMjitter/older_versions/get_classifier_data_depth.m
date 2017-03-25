function output = get_classifier_data_depth(raster,METRIC,binsize)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each

%%%%%%%%%%%%%%%%%%
iterations = 20;
%%%%%%%%%%%%%%%%%%


% Get unique depths and jitter vectors
[depths,~,ndpth] = unique([raster.AMdepth]);
[fIDs,~,nfid] = unique([raster.jitter]);


if numel(depths)<3  %skip data that can't yield a depth function
    output = 'not enough depths for neurometric data';
    return
end


% Get minimum stimulus duration
mdur = min([raster.stimDur]);


% Set up dprime matrix
dp = nan( numel(fIDs), numel(depths) );


%% Get nogo stimulus info

NGidx = ndpth==find(depths==0);
NOGO = raster(NGidx);
if numel(NOGO)>1
    NOGO = collapse_blocks(NOGO);
elseif numel(NOGO)<1
    keyboard
end

% Get nogo spiketimes, starting when sound begins
nogo_x = NOGO.x(NOGO.x>0);
nogo_y = NOGO.y(NOGO.x>0);


%% Step through the go stimuli to compare to nogo

for fid = 1:max(nfid)
    
    nYes = nan( numel(depths), 1 );
    nTrs = nan( numel(depths), 1 );
    
    for id = 1:max(ndpth)
        
        % Skip depth of 0 (nogo) and fidXdepth conditions that
        % weren't presented.
        if id==1 && depths(id)==0
            continue
        elseif id~=1 && depths(id)==0
            keyboard
        end
        if numel(raster((nfid==fid)&(ndpth==id)))<1, continue, end
        
        
        % Get all identical go trials of this type
        GO = collapse_blocks(raster((nfid==fid)&(ndpth==id)));
        
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
        [dp(fid,id),pHit(id),pFA(id)] = run_classifier(rsGO,rsNOGO,iterations);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % Get data for PY matrix fo rthis GO stim
        [nYes(id,1),nTrs(id,1)] = calculate_response_data_phys(0,pHit(id),pFA,size(rsGO,1),size(rsNOGO,1));
        
        
    end %id
    
    
    % Now get data for PY matrix for NOGO stim
    [nYes(1,1),nTrs(1,1)] = calculate_response_data_phys(1,nan,pFA,nan,size(rsNOGO,1));
    
    
    % Convert depth % to dB re 100
    depth_dB = convert_depth_proptodB(depths);
    
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Put data into output struct
    output.PYdata{fid}   = [ depth_dB'  nYes  nTrs ];
    output.dprime{fid}   = [ depth_dB(2:end)' dp(fid,2:end)'];
    output.stim{fid,1}   = fIDs(fid);
    output.stim{fid,2}   = depth_dB';
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
end %fid




end