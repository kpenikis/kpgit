function [dp,pHit,pFA] = run_classifier_Template(go_mat,ng_mat,Iterations,TemplateSize)
%---Inputs---%
%go & nogo = cell matrices of spike times for selected go and nogo trials%
%Iterations = number of times you want to run the classifier (e.g., 1000x)%
%
%  Make it so that input TemplateSize can be an integer or a fraction,
%  which determines how the template is constructed. 
% 

%---Get number of go and nogo trials---%
mxGT		=	size(go_mat,1);
mxNT		=	size(ng_mat,1);

% Set tr limit
limT = min(mxGT,mxNT);


% % % FAKE data
% % sp = randi(size(ng_mat,2),[1 floor(size(ng_mat,2)/10)]);
% % go_mat = zeros(size(go_mat));
% % go_mat(:,sp) = 1;
% % ng_mat = zeros(size(ng_mat));
% % ng_mat(:,sp) = 1;

if nargin<4
    GtmplSz      =   10;
    NtmplSz      =   10;
    keyboard
elseif TemplateSize<1 && TemplateSize>0
    GtmplSz      =   ceil(limT * TemplateSize);
    NtmplSz      =   ceil(limT * TemplateSize);
    keyboard
elseif TemplateSize>=1
    GtmplSz      =   TemplateSize;
    NtmplSz      =   TemplateSize;
    keyboard
elseif TemplateSize<0
    GtmplSz      =   limT-1;
    NtmplSz      =   limT-1;
else
    warning('detected strange input')
    keyboard
end


for ii=1:Iterations
	
	%---Get GO Template---%
    ridx        =   randperm(mxGT);
    Gidx        =   ridx(1:GtmplSz); 
    Gtemplate	=   mean(go_mat(Gidx,:),1);
    	
	%---NOGO Template---%
    ridx        =   randperm(mxNT);
    Nidx		=	ridx(1:NtmplSz); 
	Ntemplate	=	mean(ng_mat(Nidx,:),1);
	
%     figure;  plot(Gtemplate,'k')
%     hold on; plot(Ntemplate,'r')
    
	Gmn			=	nan(limT-GtmplSz,2);
	Nmn			=	nan(limT-NtmplSz,2);
	
	%---Don't grab GO template spike train vector---%
    gT          =   randperm(mxGT); %1:mxGT;
    gT          =   gT(~ismember(gT,Gidx));
    
	%---Compare GO Trials with Templates---%
	for j=1:limT-GtmplSz
		testT	=	go_mat(gT(j),:);
%         figure; plot(testT,'k')
		
		%---Compare with Templates---%
        gComp	=	pdist2(testT,Gtemplate);
        nComp	=	pdist2(testT,Ntemplate);
		Comp	=	[gComp nComp];
		Gmn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob		=	nansum(Gmn,1)/size(Gmn,1);
	
	%---Don't grab NOGO template spike train---%
    nT          =   randperm(mxNT); %1:mxNT;
    nT          =   nT(~ismember(nT,Nidx));
    
	%---Compare NOGO Trials with Templates---%
	for j=1:limT-NtmplSz
		testT	=	ng_mat(nT(j),:);
% 		figure; plot(testT,'r')
        
		%---Compare with Templates---%
        gComp	=	pdist2(testT,Gtemplate);
        nComp	=	pdist2(testT,Ntemplate);
		Comp	=	[gComp nComp];
		Nmn(j,:)	=	Comp == nanmin(Comp);
	end
	Prob2		=	nansum(Nmn,1)/size(Nmn,1);
		
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

% Correct 1s and 0s
if pHit==1
    pHit = 0.99;
elseif pHit==0
    pHit = 0.01;
end
if pFA==0
    pFA = 0.01;
elseif pFA==1
    pFA = 0.99;
end

dp			=	calculate_dprime(pHit,pFA);
% dp			=	abs(dp);

% If dp=nan and neither raster has 0 spikes; may need more iterations?
% if isnan(dp) && ~any(sum(sum(ng_mat))==0 || sum(sum(go_mat))==0)
%     keyboard
% end
%---Note: you can calculate dprime on each iteration instead and then
%average those values across all iterations. I've tried that and it's
%more/less the same results.---%
