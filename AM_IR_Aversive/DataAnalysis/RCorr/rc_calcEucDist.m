function [blockAssignment,minE] = rc_calcEucDist(T,S)
% [blockAssignment,maxR] = rc_calcR(T,S)
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
% 
%   Now calculates Euclidean Distance.


global stids

%Preallocate matrix
Rmat = nan(size(T,2),2);
Emat = nan(size(T,2),2);

for j = 1:size(T,2)
    
    % R is simply the inner dot product of S and T, divided by the
    %product of their vector norms
    if isempty(T{1,j})
        R  = nan;
        Ed = nan;
    else
%         R  = (S*T{1,j}')/(norm(S)*norm(T{1,j}));
        Ed = pdist2(S,T{1,j});
        if numel(Ed)>1
%         Ed = median(median(Ed,'omitnan'),'omitnan');
            keyboard
        end
    end
    
    %Identical calculation, but slower.
    %R = dot(S,T{1,j})/(norm(S)*norm(T{1,j}));
    
%     Rmat(j,:) = [R,j];
    Emat(j,:) = [Ed,j];
end


%Assign spike train to template level that gave highest R
% [maxR,iR] = max(Rmat(:,1));
[minE,iE] = min(Emat(:,1));

% plot(iR,iE,'ok')

%If the [min Euclidean Distance] isn't a NaN
if ~isnan(minE)
    
    %Find the index value associated with maxR
    minEindex = find(Emat(:,1) == minE);
    
    %If there's more than one index value...
    if size(minEindex,1) > 1
        
        %Randomly select one of the index 
        %values to make the assignment
        r = randi(size(minEindex,1),1);
        minEindex = minEindex(r);
    end
    
%If the max R is a NaN
elseif isnan(minE)
    
    %Make a random assignment (dont allow to choose empty stim)
    minEindex = randi(size(T,2),1);
    
end


blockAssignment = Emat(minEindex,2);


end