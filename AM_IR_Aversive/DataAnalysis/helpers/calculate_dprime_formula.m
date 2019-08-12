function dprime_mat = calculate_dprime_formula( FRmat )
% dprime_mat = calculate_dprime_formula( FR_mat )
% called by assess_this_unit
% 
% Inputs
%   FRmat (created by format_FRmat)
%   Nx4 matrix: [ stim_id  mean  std  sem ]
% 

nogo_stim = FRmat(1,1);

NOGOmean = FRmat(FRmat(:,1) == nogo_stim,2);
GOmeans = FRmat(FRmat(:,1) ~= nogo_stim,2);

NOGOstd = FRmat(FRmat(:,1) == nogo_stim,3);
GOstds = FRmat(FRmat(:,1) ~= nogo_stim,3);

%Common std
commonstds = (NOGOstd + GOstds)/2;

%Calculate dprimes
dprimes = (GOmeans-NOGOmean)./commonstds;

dprime_mat = [ FRmat(FRmat(:,1)~=nogo_stim,1) dprimes ];



end
