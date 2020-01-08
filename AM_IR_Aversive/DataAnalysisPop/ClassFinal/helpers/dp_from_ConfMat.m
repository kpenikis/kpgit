function dp = dp_from_ConfMat(ConfMat,correction,Stim)

% Do each stimulus
if nargin<3
    Stim = 1:size(ConfMat,1);
end

ConfMat(ConfMat==0) = correction;
ConfMat(ConfMat==1) = 1-correction;

dp = nan(1,numel(Stim));
for ii = 1:numel(Stim)
    ist = Stim(ii);
    dp(ii) = norminv(ConfMat(ist,ist),0,1) - norminv(mean(ConfMat(setdiff(1:end,ist),ist)),0,1);
end

end