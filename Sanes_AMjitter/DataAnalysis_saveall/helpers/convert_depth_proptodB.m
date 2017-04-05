function depth_dB = convert_depth_proptodB(depths)

% Convert depth % to dB re 100
depths(depths==0) = 0.05; depths(depths==1) = 0.95;
depth_dB = 20 .* log10(depths);

end