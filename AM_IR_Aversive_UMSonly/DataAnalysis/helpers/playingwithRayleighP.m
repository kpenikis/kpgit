% simple little script to visualize relationship between Rayleigh stat and
% its p values 

figure;
hold on

for NSpk = 0:1:300
    
    RS=[]; P=[];
    
    VSs = 0.1:0.001:0.2;
    for iVS = 1:numel(VSs)
        
        VS  = VSs(iVS);
        
        RS(iVS)     =	2*NSpk*(VS^2);         
        P(iVS)		=	RayleighP(VS,NSpk);
        
    end
    
    plot(P,RS,'.','MarkerSize',5)
    
end

xlabel('P val')
ylabel('Rayleigh')

xlim([0.000195 0.000205])

% when p = 0.0002, Rayleigh = 16.9 

