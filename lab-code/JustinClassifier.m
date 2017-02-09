function dp = RunTemplateMatching(go,nogo,Iterations)
%---Inputs---%
%go & nogo = cell matrices of spike times for selected go and nogo trials%
%Iterations = number of times you want to run the classifier (e.g., 1000x)%

%---Get number of go and nogo trials---%
mxgT		=	length(go);
mxnT		=	length(nogo);

for ii=1:Iterations
	
	%---Get GO Template---%
	Gidx		=	randi(mxgT);		%Randomly select a go trial for Go template%
	gel			=	go(:,2) == Gidx;	%go(:,2) = trial number%
	Gtemplate	=	go(gel);	
	sel			=	Gtemplate > 0;
	GtempSR		=	sum(sel);			%Go template spike count%
	
	%---NOGO Template---%
	Nidx		=	randi(mxnT);	%Randomly select a go trial for Nogo template%
	nel			=	nogoT == Nidx;	
	Ntemplate	=	nogo(nel);
	zel			=	Ntemplate > 0;
	NtempSR		=	sum(zel);			%Nogo template spike count%
	
	mn			=	nan(mxgT,2);
	Nn			=	nan(mxgT,2);
	
	%---Don't grab GO template spike train vector---%
	tel			=	go(:,2) == Gidx;
	ugTT		=	ugT(~tel);
	%---Compare GO Trials with Templates---%
	for j=1:mxgT-1
		ssel	=	goT == ugTT(j);
		Spks	=	go(ssel);
		testSR	=	Spks > 0;
		testSR	=	sum(testSR);
		
		%---Compare with Templates---%
		gComp	=	abs(testSR - GtempSR);
		nComp	=	abs(testSR - NtempSR);
		Comp	=	[gComp nComp];
		mn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob		=	nansum(mn)/length(mn);
	
	%---Don't grab NOGO template spike train---%
	ttel		=	unT == Nidx;
	unTT		=	unT(~ttel);
	%---Compare NOGO Trials with Templates---%
	for j=1:mxnT-1
		ssel	=	nogoT == unTT(j);
		Spks	=	nogo(ssel);
		testSR	=	Spks > 0;
		testSR	=	sum(testSR);
		
		%---Compare with Templates---%
		gComp	=	abs(testSR - GtempSR);
		nComp	=	abs(testSR - NtempSR);
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