function [VS,RS,P] = vectorstrength(Spktime,period,ntr)
% Vector Strength is a measure of the degree of phase-locking or
% synchronization. 
% Input Spktime must be indices, not a sparse matrix/vector.

NSpk	=	length(Spktime); %/ntr*1000/period
                                            %need this be normalized across 
                                            %different period lengths? 
                                            %ie. the average number of spikes for a 1 s sample

VS1			=	(sqrt((sum(sin(2*pi*(Spktime/period))))^2 + (sum(cos(2*pi*(Spktime/period))))^2)) / NSpk;

Hz		=	1000/period;

S		=	sum(sin(2*pi*(Spktime/1000*Hz)));
C		=	sum(cos(2*pi*(Spktime/1000*Hz)));
VS2		=	sqrt(C^2 + S^2) / length(Spktime);

if( isnan(VS1) || isnan(VS2) )
	VS	=	NaN;
else
% 	if( round(VS1,3,'decimals') ~= round(VS2,3,'decimals') )
	if( abs(VS1-VS2)>0.00001 )
		keyboard
        VS	=	VS1;
	else
		VS	=	VS1;
	end
end

RS		=	2*NSpk*(VS^2);                     % Rayleigh Statistic
P		=	RayleighP(VS,NSpk);

end
