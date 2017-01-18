function R= RayleighCrit(N,P,Step)
% RayleighCrit -- Find the critical value of the Rayleigh statistic
% given an alpha of P and N samples. Default P is .05; Step is the step 
% size that determines the precision of R. Default is .01.
% (JCM, 4/11/02)
%
% USAGE R= rayleighcrit(N,P,Step);

if nargin<3
    Step= .01;
end
if nargin<2
    P= .05;
end

for R= Step:Step:1,
    if RayleighP(R,N)<P,
        break
    end
end

return
