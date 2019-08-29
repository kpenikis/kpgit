function [blockAssignment,maxR] = rc_calcR(T,S)
% [blockAssignment,maxR] = calcR(T,S)
%   inputs
%   T: cell array of template trial data
%          (1,nblock) convolved spiketrain
%          (2,nblock) nt
%   S: convolved spiketrain for test trial
%   outputs
%   blockAssignment: best matched block from templates
%   maxR: Rcorr value of best matched block
%   
%   Calculates the normalized cross correlation (R) between a spike train 
%   and a set of template spike trains. Assigns the spike train to the 
%   template match that yields the highest R values. See (Schreiber, 2003)
%   and (Wang et al., 2007).
%
% KP, adapted from code by MLC.


global stids

%Preallocate matrix
Rmat = nan(size(T,2),2);

for j = 1:size(T,2)
    
    % R is simply the inner dot product of S and T, divided by the
    %product of their vector norms
    if isempty(T{1,j})
        R = nan;
    else
        R = (S*T{1,j}')/(norm(S)*norm(T{1,j}));
    end
    
    %Identical calculation, but slower.
    %R = dot(S,T{1,j})/(norm(S)*norm(T{1,j}));
    
    Rmat(j,:) = [R,j];
end


%Assign spike train to template level that gave highest R
maxR = max(Rmat(:,1));

%If the max R isn't a NaN
if ~isnan(maxR)
    
    %Find the index value associated with maxR
    maxRindex = find(Rmat(:,1) == maxR);
    
    %If there's more than one index value...
    if size(maxRindex,1) > 1
        
        %Randomly select one of the index 
        %values to make the assignment
        r = randi(size(maxRindex,1),1);
        maxRindex = maxRindex(r);
    end
    
%If the max R is a NaN
elseif isnan(maxR)
    
    %Make a random assignment (dont allow to choose empty stim)
    maxRindex = randi(length(stids),1);
    maxRindex = stids(maxRindex);
    
end


blockAssignment = Rmat(maxRindex,2);


end