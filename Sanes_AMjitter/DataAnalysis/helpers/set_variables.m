function [NGval,indValues,indIdx,indLabels,condIdx,condLabels] = ...
    set_variables(raster,indVar)


% Get unique depths and jitter vectors
if strcmpi(indVar,'jitter')
    raster([raster.AMdepth]==0)=[];
end
[depths,~,ridx_dpth] = unique([raster.AMdepth]);
[jitters,IA,ridx_jit] = unique([raster.jitter],'stable');

% Sort jitters to always be in increasing order
[~, isort] = sort(double(strtok(jitters,'_')));
jitters = jitters(isort);
ridx_jit_sort = nan(size(ridx_jit));
for is = 1:numel(isort)
    ridx_jit_sort(ridx_jit==IA(is),1) = isort(is);
end

switch indVar
    case 'depth'
        % Set nogo stimulus value
        NGval      = convert_depth_proptodB(0);
        % Set independent variable values and labels
        indValues  = convert_depth_proptodB(depths); % Convert depth % to dB re 100
        indIdx     = ridx_dpth;
        for il = 1:numel(depths)
            indLabels{il} = [num2str(depths(il)) '_depth'];
        end
        % Set condition variable values and labels
        condIdx    = ridx_jit_sort;
        for il = 1:numel(jitters)
            condLabels{il} = char(jitters{il});
        end
        
    case 'jitter'
        % Set nogo stimulus value
        NGval      = 0;
        % Set independent variable values and labels
        indValues  = double(strtok(jitters,'_'));
%         [indVals, isort] = sort(indValues)
        indIdx     = ridx_jit_sort;
        for il = 1:numel(indValues)
            indLabels{il} = char(jitters{il});
        end
        % Set condition variable values and labels
        condIdx    = ridx_dpth;
        for il = 1:numel(depths)
            condLabels{il} = [num2str(depths(il)) '_depth'];
        end
end

end