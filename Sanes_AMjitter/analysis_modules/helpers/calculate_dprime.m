
function dprime = calculate_dprime(pHit,pFA)

zHit	=	sqrt(2)*erfinv(2*pHit-1);
zFA		=	sqrt(2)*erfinv(2*pFA-1);

%Different way to calculate zscores
% zHit	=	norminv(pHit,0,1);
% zFA	=	norminv(pFA,0,1);

%-- Calculate d-prime
dprime = zHit - zFA ;

end