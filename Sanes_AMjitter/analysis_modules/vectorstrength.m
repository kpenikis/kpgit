function [VS,RS,P] = vectorstrength(Spks,period)
% Vector Strength is a measure of the degree of phase-locking or
% synchronization. 

NSpk	=	length(Spks);
Spktime	=	Spks;

VS1			=	(sqrt((sum(sin(2*pi*(Spktime/period))))^2 + (sum(cos(2*pi*(Spktime/period))))^2))/NSpk;

% S		=	sum(sin(2*pi*(Spktime/ISI)));
% C		=	sum(cos(2*pi*(Spktime/ISI)));

Hz		=	1000/period;

S		=	sum(sin(2*pi*(Spktime/1000*Hz)));
C		=	sum(cos(2*pi*(Spktime/1000*Hz)));
VS2		=	sqrt(C^2 + S^2) / NSpk;

if( isnan(VS1) || isnan(VS2) )
	VS	=	NaN;
else
	if( round(VS1,3,'decimals') ~= round(VS2,3,'decimals') )
		keyboard
	else
		VS	=	VS1;
	end
end

RS		=	2*NSpk*(VS^2);                     % Rayleigh Statistic
P		=	RayleighP(VS,NSpk);

end
