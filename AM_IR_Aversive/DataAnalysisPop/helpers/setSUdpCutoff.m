function q = setSUdpCutoff(Res,cutoff_val)

dp_1  =[];
dp_2  =[];
dp_3  =[];
dp_4  =[];
dp_5  =[];
% dp_all=[];

for irate = 1:size(Res,2)
    for iUn = 1:size(Res,1)
        
        if ~isempty(Res(iUn,irate).L1o)
            eval(sprintf(' dp_%i = [dp_%i; Res(iUn,irate).L1o.dprime(:,2)]; ',irate,irate));
%             dp_all = [dp_all; Res(iUn,irate).L1o.dprime(:,2)];
        end
        
    end
end


q = nan(1,5);

q(1) = quantile(dp_1,cutoff_val);
% sum(dp_1>q(1))

q(2) = quantile(dp_2,cutoff_val);
% sum(dp_2>q(2))

q(3) = quantile(dp_3,cutoff_val);
% sum(dp_3>q(3))

q(4) = quantile(dp_4,cutoff_val);
% sum(dp_4>q(4))

q(5) = quantile(dp_5,cutoff_val);
% sum(dp_5>q(5))

% q_all = quantile(dp_all,0.95);
% sum(dp_all>q_all)


end


