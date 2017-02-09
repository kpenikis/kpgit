function dprime_mat = calculatedprime(mat)
%Calculates dprime values from physiology data.
%
%Input variable mat contains either firing rate or power data, arranged as
%follows:
%
%mat = [stimulus, ave, std, sem]
%
%ML Caras Dec 2015


%If stim values are in log
if any(mat(:,1) < 1)
    nogo_target = make_stim_log(0);
else
    nogo_target = 0;
end

NOGOmean = mat(mat(:,1) == nogo_target,2);
GOmeans = mat(mat(:,1) ~= nogo_target,2);

NOGOstd = mat(mat(:,1) == nogo_target,3);
GOstds = mat(mat(:,1) ~= nogo_target,3);

%Common std
commonstds = (NOGOstd + GOstds)/2;

%Calculate dprimes
dprimes = (GOmeans-NOGOmean)./commonstds;
dprime_mat = [mat(mat(:,1) ~=nogo_target,1),dprimes];


end