function dprime_mat = calculate_dprime_formula( FR_mat )
% dprime_mat = calculate_dprime_formula( FR_mat )
%
%   This function called by format_formula_data.
% 
% Inputs
%   FR_mat (created by format_FRmat)
%   cell array, size of jitters
%   contained in each entry is an Nx4 matrix:
%   [ stim_value  mean  std  sem ]
% 

for ic = 1:numel(FR_mat)
    
    % Get data matric for this jitter 
    mat = FR_mat{ic};
    
    nogo_stim = mat(1,1);
    
    NOGOmean = mat(mat(:,1) == nogo_stim,2);
    GOmeans = mat(mat(:,1) ~= nogo_stim,2);
    
    NOGOstd = mat(mat(:,1) == nogo_stim,3);
    GOstds = mat(mat(:,1) ~= nogo_stim,3);
    
    %Common std
    commonstds = (NOGOstd + GOstds)/2;
    
    %Calculate dprimes
    dprimes = (GOmeans-NOGOmean)./commonstds;
    
    dprime_mat{ic} = [ mat(mat(:,1)~=nogo_stim,1) dprimes ];
    
end

end
