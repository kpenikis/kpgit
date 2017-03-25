function P= RayleighP(R,N)
% RayleighP -- Find the probability of vector strength R given N measurements.
% JCM, from Mardia
%
% USAGE P= rayleighp(R,N);

K= N*R^2;

P= exp(-K)*(1+(2*K-K^2)/(4*N) - (24*K-132*K^2+76*K^3 - 9*K^4)/(288*N^2));
