function [PYdata_out,dprime_out,stim] = get_classifier_data(raster,METRIC,binsize,indVar)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
%
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each

global subject

%%%%%%%%%%%%%%%%%%
iterations = 1000;
%%%%%%%%%%%%%%%%%%


% Set independent and condition variables
[NGval,indValues,indIdx,indLabels,condIdx,condLabels] = set_variables(raster,indVar);

if numel(indValues)<3  %skip data that can't yield a depth function -->should not be called though
    keyboard
    return
end

% Get minimum stimulus duration
mdur = min([raster.stimDur]);

rng('shuffle')

%% Get nogo stimulus info

NGidx = indIdx==find(indValues==NGval);
NOGO = raster(NGidx);
if numel(NOGO)>1
    if strcmp(indVar,'jitter') && any([NOGO.AMdepth]==0)
        NOGO([NOGO.AMdepth]==0)=[]; % Get rid of the unmodulated stim for now. 
                                    % Eventually may want to discriminate
                                    % against it too.
    else
        keyboard
%         NOGO = collapse_blocks(NOGO);
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
        
        if numel(raster((condIdx==ic)&(indIdx==ii)))<1, continue,
        elseif numel(raster((condIdx==ic)&(indIdx==ii)))>1, keyboard, end
        
        
        % Get all identical go trials of this type
        GO = raster((condIdx==ic)&(indIdx==ii));
        
        % Get spiketimes for go stimulus, starting when sound begins
        go_x = GO.x(GO.x>0);
        go_y = GO.y(GO.x>0);
        
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
        
        % Get trial vectors to discriminate
        if strcmp(METRIC,'VS')
            
            for iy = 1:max(go_y)
                trials = randperm(max(go_y));
                rsGO(iy,1) = calc_VSRS(GO,subject,trials(1:10));
            end
            for iy = 1:max(go_y)
                trials = randperm(max(go_y));
                rsNOGO(iy,1) = calc_VSRS(NOGO,subject,trials(1:10));
            end
            
        elseif strcmp(METRIC,'Corr')
            
            % For the GO stimulus
            [Rs,~,sh,~,~,~]    = corr_spks(GO,subject);
            [~,peakLag]=findpeaks(Rs,'MinPeakProminence',0.005);
            if numel(peakLag)<1
                [~,peakLag]=findpeaks(Rs,'MinPeakProminence',0.003);
            end
            if numel(peakLag)<1
                disp('weak correlations. setting lag to 0')
                peakLag = find(sh==0);
            elseif numel(peakLag)>1
                [~,leastshift] = min(abs(peakLag-find(sh==0)));
                peakLag = peakLag(leastshift);
            end
            [~,~,~,rsGO,~,~]   = corr_spks(GO,subject,sh(peakLag));
            
            % For the NOGO stimulus
            [Rs,~,sh,~,~,~]    = corr_spks(NOGO,subject);
            peakLag=[];
            [~,peakLag]=findpeaks(Rs,'MinPeakProminence',0.005);
            if numel(peakLag)<1
                [~,peakLag]=findpeaks(Rs,'MinPeakProminence',0.003);
            end
            if numel(peakLag)<1
                disp('weak correlations. setting lag to 0')
                peakLag = find(sh==0);
            elseif numel(peakLag)>1
                [~,leastshift] = min(abs(peakLag-find(sh==0)));
                peakLag = peakLag(leastshift);
            end
            [~,~,~,rsNOGO,~,~] = corr_spks(NOGO,subject,sh(peakLag));
            
        else
            
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
    
    % Correct infinite dprime vals
    for inf_idx = find(dprime==Inf)
        if pHit(inf_idx)==1
            pHit(inf_idx) = 0.999;
        end
        if pFA(inf_idx)==0
            pFA(inf_idx) = 0.001;
        end
        corrected_dp = calculate_dprime(pHit(inf_idx),pFA(inf_idx));
        dprime(inf_idx) = corrected_dp;
    end
    
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Put data into output struct
    PYdata_out{ic}   = [ indValues'  nYes  nTrs ];
    dprime_out{ic}   = [ indValues(2:end)' dprime(ic,2:end)'];
    stim{ic,1}   = condLabels{ic};
    stim{ic,2}   = indLabels';
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
end %ic




end