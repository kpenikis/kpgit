function [blocks,group_blocks] = set_grouped_blocks(blocks)
% Groups blocks from identical recordings that were interrupted and should
% be combined manually.


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~
group_blocks = [90 89];  %only 2 at a time for now
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for ic = size(group_blocks,1)
    blocks(blocks==group_blocks(ic,2)) = [];
end



end