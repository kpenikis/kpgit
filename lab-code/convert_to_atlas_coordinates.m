
function [AP_new, DV_new, MLcoord] = convert_to_atlas_coordinates(AP,ML,SD)

% Convert my coordinates to slightly smaller brain
% AP = AP*(4.45/5);
% ML = ML*(4.45/5);
% SD = SD*(4.45/5);

% Convert surgical depth to DV axis depth
DV = SD * cos(deg2rad(25));
MLcoord = SD * sin(deg2rad(25)) + ML;

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






