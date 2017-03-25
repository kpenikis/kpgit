function IRseq_USE = make_IR_segments

AMrates = 2.^[1:6];

indices = ...
   [1 2 3 4 5 6;...
    6 5 4 3 2 1;...
    2 4 6 1 3 5;...
    5 3 1 6 4 2;...
    3 6 2 5 1 4;...
    4 1 5 2 6 3];

for is = 1:100
    seed1 = randperm(6);
    IRseqs(:,:,is) = seed1(indices);
end

stepsize2 = diff(IRseqs,2,2)~=0;

movement_metric=[];
for is = 1:size(IRseqs,3)
movement_metric(is) = sum(sum(stepsize2(:,:,is),1),2);
end

max_mvt = max(movement_metric);
max_idxs=[];
max_idxs = find(movement_metric==max_mvt);

IRseqs = IRseqs(:,:,max_idxs);

for is = 1:size(IRseqs,3)
    this_set = IRseqs(:,:,is);
    this_set = this_set(this_set(:,1)~=3 & this_set(:,1)~=4,:);

    [i,j]=find(this_set==3);
    ij = sub2ind(size(this_set),i,j-1);
    pre3 = this_set(ij);
    
    [i,j]=find(this_set==4);
    ij = sub2ind(size(this_set),i,j-1);
    pre4 = this_set(ij);
    
    if sum(pre3-3)==0 && sum(pre4-4)==0
        IRseq_USE = this_set;
        break
    end
    
end

for is = 1:4
    IRsegmentA = AMrates(IRseq_USE(1,:));
    IRsegmentB = AMrates(IRseq_USE(3,:));
    IRsegmentC = AMrates(IRseq_USE(2,:));
    IRsegmentD = AMrates(IRseq_USE(4,:));
end

end


