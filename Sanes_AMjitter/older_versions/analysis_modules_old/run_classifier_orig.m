function dp = run_classifier_orig(go_x,go_y,nogo_x,nogo_y,Iterations)
%---Inputs---%
%go & nogo = cell matrices of spike times for selected go and nogo trials%
%Iterations = number of times you want to run the classifier (e.g., 1000x)%

%---Get number of go and nogo trials---%
mxgT		=	max(go_y);
mxnT		=	max(nogo_y);

for ii=1:Iterations
	
	%---Get GO Template---%
	Gidx		=	randi(mxgT);		%Randomly select a go trial for Go template%
	Gtemplate	=	go_x(go_y==Gidx);	
	GtempSR		=	length(Gtemplate);	%Go template spike count%
	
	%---NOGO Template---%
	Nidx		=	randi(mxnT);	    %Randomly select a go trial for Nogo template%
	Ntemplate	=	nogo_x(nogo_y==Nidx);
	NtempSR		=	length(Ntemplate);  %Nogo template spike count%
	
	mn			=	nan(mxgT-1,2);
	Nn			=	nan(mxnT-1,2);
	
	%---Don't grab GO template spike train vector---%
    gT          =   1:mxgT;
    gT          =   gT(gT~=Gidx);
	%---Compare GO Trials with Templates---%
	for j=1:mxgT-1
		Spks	=	go_x(go_y==gT(j));
        testSR	=	length(Spks);
		
		%---Compare with Templates---%
% 		gComp_a	=	abs(testSR - GtempSR);
        gComp	=	pdist2(testSR,GtempSR);
%         nComp_a	=	abs(testSR - NtempSR);
        nComp	=	pdist2(testSR,NtempSR);
		Comp	=	[gComp nComp];
		mn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob		=	nansum(mn)/length(mn);
	
	%---Don't grab NOGO template spike train---%
    nT          =   1:mxnT;
    nT          =   nT(nT~=Nidx);
	%---Compare NOGO Trials with Templates---%
	for j=1:mxnT-1
		Spks	=	nogo_x(nogo_y==nT(j));
		testSR	=	length(Spks);
		
		%---Compare with Templates---%
% 		gComp_a	=	abs(testSR - GtempSR);
        gComp	=	pdist2(testSR,GtempSR);
% 		nComp_a	=	abs(testSR - NtempSR);
        nComp	=	pdist2(testSR,NtempSR);
		Comp	=	[gComp nComp];
		Nn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob2		=	nansum(Nn)/length(Nn);
		
	if( ii == 1 )
		MN		=	Prob;
		NN		=	Prob2;
	else
		MN		=	[MN;Prob];
		NN		=	[NN;Prob2];
	end
end
p			=	nanmean(MN,1);
pp			=	nanmean(NN,1);
pHit		=	p(1);
pFA			=	pp(1);
dp			=	calculatedprime(pHit,pFA);
dp			=	abs(dp);
%---Note: you can calculate dprime on each iteration instead and then
%average those values across all iterations. I've tried that and it's
%more/less the same results.---%

function dprime = calculatedprime(pHit,pFA)
zHit	=	sqrt(2)*erfinv(2*pHit-1);
zFA		=	sqrt(2)*erfinv(2*pFA-1);
%Differe way to calculate dprime%
% zHit	=	norminv(pHit,0,1);
% zFA		=	norminv(pFA,0,1);
%-- Calculate d-prime
dprime = zHit - zFA ;