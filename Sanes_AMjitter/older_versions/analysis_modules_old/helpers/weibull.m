function p = weibull(a,B,x,g)
p = 1 - g*( exp( -(x/a).^B )) ;
end