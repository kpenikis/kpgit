function [FRmat,stim,ntrials] = format_FRmat( raster, respdur, indVar )
% called by format_formula_data


% Set independent and condition variables
[NGval,indValues,indIdx,indLabels,condIdx,condLabels] = set_variables(raster,indVar);

if numel(indValues)<3  %skip data that can't yield a nm function
    FRmat = 'not enough depths for neurometric data';
    return
end


%% Get NOGO first (because there is one for all jitters)
data = raster(indIdx(indValues==NGval));

if ischar(respdur) %use full stimulus duration
    dur = data.stimDur; 
else               %use specified duration from onset
    dur = respdur;
end
y = data.y(data.x>0 & data.x<dur);
x = data.x(data.x>0 & data.x<dur);

fr_vec = zeros(max(y),1);
for iy = 1:max(y)
    fr_vec(iy) = numel(x(y==iy))/(dur/1000);
end
ng_avgFR = mean(fr_vec); %sp/s
ng_stdFR = std(fr_vec);
ng_semFR = ng_stdFR/sqrt(max(y) - 1);
ng_ntr   = max(y);


%% Now go through GO stimuli
for ic = 1:max(condIdx)
    
    % Create stim output
    stim{ic,1}   = condLabels{ic};
    stim{ic,2}   = indLabels';
    
    %FRmat = [stimulus, aveFR, stdFR, semFR]
    FRmat{ic} = nan( numel(indValues), 4 );
    % Fill in NOGO data
    FRmat{ic}(1,:) = [indValues(1) ng_avgFR ng_stdFR ng_semFR];
    ntrials{ic}(1,1) = ng_ntr;
    
    for ii = 2:max(indIdx)
        
        data = raster((condIdx==ic)&(indIdx==ii));
        if isempty(data), continue, end
        
        x=[]; y = [];
        if ischar(respdur) %use full stimulus duration
            dur = data.stimDur;
        else               %use specified duration from onset
            dur = respdur;
        end
        y = data.y(data.x>0 & data.x<dur);
        x = data.x(data.x>0 & data.x<dur);
        
        
        fr_vec = zeros(max(y),1);
        for iy = 1:max(y)
            fr_vec(iy) = numel(x(y==iy))/(dur/1000);
        end
        avgFR = mean(fr_vec); %sp/s
        stdFR = std(fr_vec);
        semFR = ng_stdFR/sqrt(max(y) - 1);
        
        % Fill in GO data
        FRmat{ic}(ii,:) = [indValues(ii) avgFR stdFR semFR];
        ntrials{ic}(ii,1) = max(y);
        
    end
end




end