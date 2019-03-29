function [rn,rt] = noisecorr(N1,N2)
%---Average FR per stimulus---%
FR1             =   cell2mat(N1.SpkRate);
FR2             =   cell2mat(N2.SpkRate);
%---FR per stimulus trial---%
SpkVec1         =   N1.SpkVec;
SpkVec2         =   N2.SpkVec;
%---R-noise---%
rt              =   corr(FR1,FR2);

Nrates          =   length(FR1);
rn              =   nan(1,Nrates);
for i=1:Nrates
    
    n1          =   SpkVec1{i};
    n2          =   SpkVec2{i};
    
    rn(1,i)     =   corr(n1,n2);
    
end