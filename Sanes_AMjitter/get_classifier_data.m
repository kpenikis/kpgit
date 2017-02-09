function output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
% output = get_classifier_data(raster,METRIC,binsize,merge_blocks)
% 
%  First organizes all unique stimuli, then steps through to get AM depth
% detection discriminability data for each 

%%%%%%%%%%%%%%%%%%
iterations = 1000;
%%%%%%%%%%%%%%%%%%


if ~merge_blocks
    % Get blocks for loop and designate which ones to combine
    [blocks,group_blocks] = set_grouped_blocks(unique([raster.block]));
else
    % this case if you do want to merge across blocks (i.e. depth
    % discrimination when blocks of 75% depth and 41% depth are presented
    % separately).
    keyboard
    blocks = 1;
end

% Make separate figures for each param except jitter and depth
for ibk = blocks
    
    if ~merge_blocks
        % Get data for these blocks
        bk_raster = raster([raster.block]==ibk);
        
        % Add data from blocks set to combine, if needed
        if ~isempty(group_blocks(:,1)==ibk)
            bk_raster = [bk_raster raster([raster.block]== group_blocks(group_blocks(:,1)==ibk,2) )];
            bk_str = [num2str(ibk) num2str(group_blocks(group_blocks(:,1)==ibk,2))];
        else
            bk_str = num2str(ibk);
        end
        
    else
        bk_raster = raster;
        bk_str = 'all';
    end
    
    
    % Find unique stimuli based on other parameters
    [LP_HP_dB_rate,~,np] = unique([bk_raster.HP; bk_raster.LP; bk_raster.dB; bk_raster.AMrate]','rows');
    
    
    for ip = 1:max(np)
        
        param_raster = bk_raster(np==ip);
        
        % Separate figures by behavioral state
%         behavstates = {'D' 'P' 'A'};
        for ib = {'D' 'P' 'A'}
            
            stim = param_raster(strcmp({param_raster.behaving},ib));
            
            % Get unique depths and jitter vectors
            [depths,~,ndpth] = unique([stim.AMdepth]);
            [fIDs,~,nfid] = unique({stim.stimfn},'stable');
            
            if numel(depths)<3  %skip data that can't yield a depth function
                continue
            end
            
            % Right now, data is saved to output struct in a way that is not
            % compatible with multiple other stimulus params.
            if max(np)>1, keyboard, end
            
            mdur = min([stim.stimDur]);
            dp = nan( numel(fIDs), numel(depths) );
            
            % Get nogo stimulus info
            NGidx = ndpth==find(depths==0);
            NOGO = stim(NGidx);
            if numel(NOGO)>1
                NOGO = collapse_blocks(NOGO);
            elseif numel(NOGO)<1
                keyboard
            end
            
            % Get nogo spiketimes, starting when sound begins
            nogo_x = NOGO.x(NOGO.x>0);
            nogo_y = NOGO.y(NOGO.x>0);
            
            %%
            % Step through the go stimuli to compare to nogo
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
                    if numel(stim((nfid==fid)&(ndpth==id)))<1, continue, end
                    
                    
                    % Get all identical go trials of this type
                    GO = collapse_blocks(stim((nfid==fid)&(ndpth==id)));
                    
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
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    % Run classifier
                    [dp(fid,id),pHit(id),pFA(id)] = run_classifier(rsGO,rsNOGO,iterations);
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    
                    % Get data for PY matrix fo rthis GO stim
                    [nYes(id,1),nTrs(id,1)] = calculate_response_data_phys(0,pHit(id),pFA,length(rsGO),length(rsNOGO));
                                        
                    
                end %id
                
                
                % Now get data for PY matrix for NOGO stim
                [nYes(1,1),nTrs(1,1)] = calculate_response_data_phys(1,nan,pFA,nan,length(rsNOGO));
                
                
                % Convert depth % to dB re 100
                depths(depths==0) = 0.05; depths(depths==1) = 0.95;
                depth_dB = 20 .* log(depths);
                
                
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                % Put data into output struct
                output.PYdata{fid}   = [ depth_dB'  nYes  nTrs ];
                output.dprime{fid}   = [ depth_dB(2:end)' dp(fid,2:end)'];
                output.stim{fid,1}   = fIDs{fid};
                output.stim{fid,2}   = depth_dB';
                output.stim{fid,3}   = ib{1};
                output.stim{fid,4}   = LP_HP_dB_rate(ip,:);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
            end %fid
            
            
                        
            
            
        end
    end
end
end