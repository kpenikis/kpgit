function [VS,RS,P,mu] = vectorstrength_meanphase(Spktime,period)
% [VS,RS,P] = vectorstrength_meanphase(Spktime,period)
%   Vector Strength is a measure of the degree of phase-locking or
%   synchronization. 
%   Input Spktime must be indices, not a sparse matrix/vector.
% 

if isempty(Spktime)
    VS = 0;
    RS = 0;
    P  = 1;
    mu = nan;
    return
end


NSpk	=	length(Spktime);  %/ntr*1000/period
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
	if( abs(VS1-VS2)>0.001 )
		keyboard
        VS	=	VS1;
	else
		VS	=	VS1;
	end
end

RS		=	2*NSpk*(VS^2);                     % Rayleigh Statistic
P		=	RayleighP(VS,NSpk);


% Calculate mean phase 

if VS>0
    
    % output of atan: [ -pi/2, pi/2 ]
    % Must correct based on quadrant to get mu in range [ 0, 2*pi ]
    
    if  C<0                 % Q2 or Q3  : flip across y axis, output between [pi/2,3*pi/2]
        mu      =   atan(S/C) + pi;
    elseif  C>=0 && S>=0    % Q1  :  output between [0,pi/2]
        mu      =   atan(S/C);   
    elseif  C>=0 && S<0     % Q4  :  output between [-pi/2,0], shift to [3*pi/2,2*pi]
        mu      =   atan(S/C) + 2*pi;
    else
        keyboard
    end
else
    if C==S
        mu      =   atan(1);
    end
    keyboard
end


end
