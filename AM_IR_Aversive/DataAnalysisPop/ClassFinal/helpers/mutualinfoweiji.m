%
%  0.5  0.4  0.1 ---> add up to 1, represent p(Y|X=1)
%  0    0.7  0.3 ---> add up to 1, represent p(Y|X=2)
%  0.2  0.2  0.6 ---> add up to 1, represent p(Y|X=3)
% 
% Now let's say that 
% p(X) is the distribution of the stimuli X 
% p(Y) is the distribution of the responses Y regardless of the true stimulus, i.e.
% 
% p(Y) = sum_X p(X) p(Y|X)
% 
% 1) For each X and Y, calculate p(Y|X) log p(Y|X)
% 2) Sum over Y.
% 3) Do a weighted average over X, with weights equal to p(X)
% 4) Subtract p(Y) log p(Y) summed over Y



ConfMat

x = (1:8)';
p_x = ones(size(x))./numel(x);

p_yx = nan(1,numel(x));
for ist = x'
    p_yx(ist) = sum(ConfMat(ist,:));
end


p_y = p_x'.*log(p_yx);


sum_y = sum(ConfMat,1);
