
function [AP_new, DV_new, MLcoord] = convert_to_atlas_coordinates(AP,ML,SD)

% Convert lab coordinates to slightly smaller brain
% (tried scaling to see if this could reduce the discrepency
%  between our A1 coords and where that brings you to in the atlas)
% AP = AP*(4.45/5);
% ML = ML*(4.45/5);
% SD = SD*(4.45/5);

% Convert surgical depth to DV axis depth
electrode_angle = 25;%deg
DV = SD * cos(deg2rad(electrode_angle));
MLcoord = SD * sin(deg2rad(electrode_angle)) + ML;

% Get angle for coordinate shift
% Based on the distance between lambda and bregma, and the amount that
% bregma falls ventral to lambda when skull aligned using occipital crest.
lambda_bregma_dist = 4.45;
bregma_drop = 1.5;
shift_angle = atan(bregma_drop/lambda_bregma_dist); %rad

% Make the conversion
H = sqrt(AP^2 + DV^2);
x = pi/2 - shift_angle - acos(AP/H);

AP_new = H * sin(x);
DV_new = H * cos(x);






