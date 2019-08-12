function [dp,pHit,pFA] = run_classifier_SingleTrials(go_mat,ng_mat,Iterations)
%---Inputs---%
%go & nogo = cell matrices of spike times for selected go and nogo trials%
%Iterations = number of times you want to run the classifier (e.g., 1000x)%

%---Get number of go and nogo trials---%
mxgT		=	size(go_mat,1);
mxnT		=	size(ng_mat,1);

for ii=1:Iterations
	
	%---Get GO Template---%
	Gidx		=	randi(mxgT);		%Randomly select a go trial for Go template%
	Gtemplate	=   go_mat(Gidx,:);  
% 	GtempSR		=	length(Gtemplate);	%Go template spike count%
	
	%---NOGO Template---%
	Nidx		=	randi(mxnT);	    %Randomly select a go trial for Nogo template%
	Ntemplate	=	ng_mat(Nidx,:);
	
	mn			=	nan(mxgT-1,2);
	Nn			=	nan(mxnT-1,2);
	
	%---Don't grab GO template spike train vector---%
    gT          =   1:mxgT;
    gT          =   gT(gT~=Gidx);
	%---Compare GO Trials with Templates---%
	for j=1:mxgT-1
		testT	=	go_mat(gT(j),:);
		
		%---Compare with Templates---%
        gComp	=	pdist2(testT,Gtemplate);
        nComp	=	pdist2(testT,Ntemplate);
		Comp	=	[gComp nComp];
		mn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob		=	nansum(mn)/length(mn);
	
	%---Don't grab NOGO template spike train---%
    nT          =   1:mxnT;
    nT          =   nT(nT~=Nidx);
	%---Compare NOGO Trials with Templates---%
	for j=1:mxnT-1
		testT	=	ng_mat(nT(j),:);
		
		%---Compare with Templates---%
        gComp	=	pdist2(testT,Gtemplate);
        nComp	=	pdist2(testT,Ntemplate);
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
dp			=	calculate_dprime(pHit,pFA);
dp			=	abs(dp);
%---Note: you can calculate dprime on each iteration instead and then
%average those values across all iterations. I've tried that and it's
%more/less the same results.---%
