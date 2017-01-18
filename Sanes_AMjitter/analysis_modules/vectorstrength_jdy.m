function [VS,RS,P] = vectorstrength_jdy(Spks,ISI)

% Vector Strength is a measure of the degree of phase-locking or
% synchronization. (Computed over entire stimulus interval, including rise
% and fall time, and not accounting for latency. This should be changed
% later).

NSpk	=	length(Spks);
Spktime	=	Spks;

VS1			=	(sqrt((sum(sin(2*pi*(Spktime/ISI))))^2 + (sum(cos(2*pi*(Spktime/ISI))))^2))/NSpk;

% S		=	sum(sin(2*pi*(Spktime/ISI)));
% C		=	sum(cos(2*pi*(Spktime/ISI)));

Hz		=	1000/ISI;
Hz		=	Hz*0.001;
S		=	sum(sin(2*pi*(Spktime*Hz)));
C		=	sum(cos(2*pi*(Spktime*Hz)));
VS2		=	sqrt(C^2 + S^2) / NSpk;
if( isnan(VS1) || isnan(VS2) )
	VS	=	NaN;
else
	if( dround(VS1) ~= dround(VS2) )
		keyboard
	else
		VS	=	VS1;
	end
end

RS		=	2*NSpk*(VS^2);                     % Rayleigh Statistic

P		=	RayleighP(VS,NSpk);
